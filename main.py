import sys
import subprocess
import requests
from pathlib import Path
from tqdm import tqdm

from config import PipelineConfig, parse_args, setup_logging
from modules.ena_downloader import ENADownloader
from modules.bwa_mapper import BWAMapper
from modules.softclipper import SoftClipper
from modules.bam_processor import BAMProcessor
from modules.analyzers import MapDamageAnalyzer, QualimapAnalyzer, HaplotypeCaller


def _cleanup_intermediate_file(path: Path, logger):
    """Remove intermediate files safely while logging the result."""
    if not path:
        return
    try:
        if path.exists():
            path.unlink()
            logger.info(f"Removed intermediate file: {path}")
    except Exception as exc:  # pragma: no cover - cleanup failures shouldn't stop the pipeline
        logger.warning(f"Failed to remove intermediate file {path}: {exc}")

def main():
    # 設定とログの初期化
    args = parse_args()
    config = PipelineConfig(args)
    
    # ログファイルのパスを指定（どのファイルにログを残すか明示）
    log_file = config.logs_dir / f"pipeline_{config.project_accession}.log"
    logger = setup_logging(log_file=log_file)
    
    logger.info(f"Starting complete pipeline for project {config.project_accession}")
    
    # Check reference genome
    if not config.reference_genome.exists():
        logger.error(f"Reference genome not found: {config.reference_genome}")
        sys.exit(1)
    
    # Create reference index if needed
    if not config.reference_genome.with_suffix('.fai').exists():
        logger.info("Creating reference genome index")
        subprocess.run(["samtools", "faidx", str(config.reference_genome)], check=True)
    
    # モジュールの初期化
    ena_downloader = ENADownloader(config)
    bwa_mapper = BWAMapper(config)
    softclipper = SoftClipper(config)
    bam_processor = BAMProcessor(config)
    mapdamage_analyzer = MapDamageAnalyzer(config)
    qualimap_analyzer = QualimapAnalyzer(config)
    haplotypecaller = HaplotypeCaller(config)
    
    session = requests.Session()
    
    # Step 1: Download data from ENA
    logger.info("Step 1: Downloading data from ENA")
    response_data = ena_downloader.get_api_response(config.project_accession, session)
    sample_to_ftp_urls = ena_downloader.parse_response_data(response_data)
    
    # サンプルごとの処理をプログレスバーで表示
    for sample_acc, ftp_urls in tqdm(sample_to_ftp_urls.items(), desc="progress", unit="sample"):
        logger.info(f"Processing sample: {sample_acc}")
        
        # Download FASTQ files
        fastq_files = ena_downloader.download_sample_data(sample_acc, ftp_urls)
        
        if not fastq_files:
            logger.warning(f"No FASTQ files downloaded for {sample_acc}")
            continue
        
        # Step 2: BWA mapping
        bam_file = bwa_mapper.run_mapping_pipeline(sample_acc, fastq_files)
        if not bam_file:
            logger.error(f"Mapping failed for {sample_acc}")
            continue
        
        # Step 3: Softclipping
        softclipped_bam = softclipper.run_softclipping(sample_acc, bam_file)
        if not softclipped_bam:
            logger.error(f"Softclipping failed for {sample_acc}")
            continue
        _cleanup_intermediate_file(bam_file, logger)

        # Step 4: BAM processing (sort, dedup, index)
        dedup_bam = bam_processor.run_bam_processing(sample_acc, softclipped_bam)
        if not dedup_bam:
            logger.error(f"BAM processing failed for {sample_acc}")
            continue

        # Step 5: mapDamage
        mapdamage_result = mapdamage_analyzer.run_mapdamage(sample_acc, softclipped_bam)
        if not mapdamage_result:
            logger.error(f"MapDamage analyzing failed for {sample_acc}")
            continue

        bam_dir = dedup_bam.parent.parent / "bam_files"
        if bam_dir.exists():
            for pattern in ("*.bam", "*.bai", "*.truncated"):
                for intermediate in bam_dir.glob(pattern):
                    _cleanup_intermediate_file(intermediate, logger)

        # Step 6: Qualimap
        qualimap_result = qualimap_analyzer.run_qualimap(sample_acc, dedup_bam)
        if not qualimap_result:
            logger.error(f"Qualimap analyzing failed for {sample_acc}")
            continue
        
        # Step 7: HaplotypeCaller
        vcf_file = haplotypecaller.run_haplotypecaller(sample_acc, dedup_bam)
        if not vcf_file:
            logger.error(f"HaplotypeCaller failed for {sample_acc}")
            continue

        # Cleanup large intermediate files once downstream analyses are finished
        dedup_bam_index = Path(str(dedup_bam) + ".bai")
        _cleanup_intermediate_file(dedup_bam_index, logger)
        _cleanup_intermediate_file(dedup_bam, logger)
        _cleanup_intermediate_file(softclipped_bam, logger)

        logger.info(f"Completed processing for {sample_acc}")
    
    logger.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()