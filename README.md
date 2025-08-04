# auto-app2

## 概要

**auto-app2** は NGS（次世代シーケンス）データ解析を自動化する Python 製パイプラインです。
ENA からのデータ取得、BWA によるマッピング、Soft clipping、Picard や samtools による BAM 処理、mapDamage/Qualimap による QC、GATK HaplotypeCaller を用いた variant calling までを一括で実行します。古 DNA やゲノム解析の再現性の高いワークフロー構築を支援します。

- **対応 OS:** macOS, Linux
- **主な機能:** ENA からの FASTQ 取得、BWA マッピング、Soft clipping、BAM sort/duplicate removal、mapDamage・Qualimap・HaplotypeCaller 解析、VCF 出力

---

## 目次

- [特徴](#特徴)
- [ディレクトリ構成](#ディレクトリ構成)
- [必要条件](#必要条件)
- [セットアップ](#セットアップ)
  - [1. Conda での環境構築](#1-conda-での環境構築)
  - [2. Docker での実行](#2-docker-での実行)
  - [3. Apple Silicon の注意点](#3-apple-silicon-の注意点)
- [使い方](#使い方)
  - [コマンドライン引数](#コマンドライン引数)
  - [実行例](#実行例)
- [ワークフロー詳細](#ワークフロー詳細)
- [各種ディレクトリ・ファイルの説明](#各種ディレクトリファイルの説明)
- [トラブルシューティング・FAQ](#トラブルシューティングfaq)

---

## 特徴

- **自動化**: データ取得から解析まで一括実行
- **再現性**: コマンドライン引数でパラメータを管理
- **拡張性**: モジュール構造で機能追加が容易

---

## ディレクトリ構成

```
auto-app2/
├── config.py                # 設定と引数解析
├── environment.yml          # Conda 環境定義
├── main.py                  # パイプライン実行スクリプト
├── modules/                 # 個別モジュール
│   ├── analyzers.py         # mapDamage・Qualimap・HaplotypeCaller
│   ├── bam_processor.py     # BAM sort / dedup / index
│   ├── bwa_mapper.py        # BWA + AdapterRemoval
│   ├── ena_downloader.py    # ENA から FASTQ ダウンロード
│   ├── softclipper.py       # Soft clipping 処理
│   └── __init__.py
├── data/                    # 実行時に作成されるディレクトリ
│   ├── raw_data/            # ダウンロードした FASTQ
│   ├── reference/           # 参照ゲノム (手動で配置)
│   ├── results/             # 解析結果
│   ├── logs/                # 実行ログ
│   └── temp/                # 一時ファイル
├── requirements.txt         # 追加 Python 依存関係
└── README.md
```

---

## 必要条件

- Python 3.8 以上
- BWA、samtools、AdapterRemoval、Picard、mapDamage、Qualimap、GATK などの外部ツール
- インターネット接続 (ENA からのデータ取得に使用)
- 解析用の参照ゲノム FASTA ファイル

---

## セットアップ

### 1. Conda での環境構築

```sh
# Conda のインストール (未導入の場合)
# macOS
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install --cask miniconda

# Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# 環境作成
conda env create -f environment.yml
conda activate auto-app2-env

# 追加 Python パッケージ
pip install -r requirements.txt
```

必要に応じて AdapterRemoval・Picard・mapDamage・Qualimap・GATK などを conda または別途インストールしてください。

### 2. Docker での実行

Docker/Docker Compose を利用する場合:

```sh
docker compose build
docker compose run --rm auto-app2 python main.py --project_accession PRJEB19970 --reference_genome /app/data/reference/equCab3.fa
```

### 3. Apple Silicon の注意点

Apple Silicon (M1～M4) では一部のツールが conda からインストールできない場合があります。その際は Homebrew 等で手動インストールしてください。

```sh
brew install bwa
brew install brewsci/bio/adapterremoval
brew install samtools
```

---

## 使い方

### コマンドライン引数

`main.py` は次の引数を受け取ります。

| 引数 | 説明 | デフォルト |
| ---- | ---- | ---------- |
| `--project_accession` | ENA のプロジェクト ID | `PRJEB19970` |
| `--base_dir` | データ保存先のベースディレクトリ | `./data` |
| `--reference_genome` | 参照ゲノム FASTA ファイル | `./data/reference/equCab3.fa` |
| `--workers` | ENA ダウンロードの並列数 | `4` |
| `--threads` | 各種解析のスレッド数 | `20` |
| `--java_mem` | Java ツール用メモリ設定 | `10g` |

### 実行例

```sh
# 基本的な実行例
python main.py \
  --project_accession PRJEB19970 \
  --reference_genome ./data/reference/equCab3.fa
```

---

## ワークフロー詳細

1. **ENA から FASTQ をダウンロード** (`modules/ena_downloader.py`)
2. **BWA + AdapterRemoval でマッピング** (`modules/bwa_mapper.py`)
3. **Soft clipping の適用** (`modules/softclipper.py`)
4. **BAM sort / dedup / index** (`modules/bam_processor.py`)
5. **mapDamage によるダメージ解析** (`modules/analyzers.py`)
6. **Qualimap による品質評価** (`modules/analyzers.py`)
7. **HaplotypeCaller による Variant Calling** (`modules/analyzers.py`)

各ステップは `main.py` から自動的に呼び出され、結果は `data/results/<ProjectID>/<SampleID>/` 以下に保存されます。

---

## 各種ディレクトリ・ファイルの説明

- `config.py`: コマンドライン引数とログ設定
- `environment.yml`: Conda で使用する依存パッケージ
- `requirements.txt`: 追加の Python パッケージ
- `main.py`: パイプラインのエントリーポイント
- `modules/`: 各機能を実装した Python モジュール
- `data/raw_data/`: ダウンロードした FASTQ
- `data/reference/`: 参照ゲノム (ユーザーが配置)
- `data/results/<ProjectID>/<SampleID>/`: サンプルごとの解析結果 (BAM, VCF, QC レポートなど)
- `data/logs/`: 実行ログ
- `data/temp/`: 一時ファイル

---

## トラブルシューティング・FAQ

### Q. 依存パッケージのエラーが出る  
A. `conda activate auto-app2-env` を実行したか確認し、`pip install -r requirements.txt` を再実行してください。

### Q. 外部コマンドが見つからない (bwa, AdapterRemoval など)  
A. Apple Silicon Mac では Homebrew などで手動インストールしてください。

### Q. データがダウンロードできない  
A. プロジェクト ID やネットワーク接続を確認してください。

### Q. BAM/VCF などの出力がない  
A. `data/logs/` に出力されるログを確認し、エラー内容を特定してください。

### Q. Qualimap 実行時に `Unrecognized VM option 'MaxPermSize=1024m'` と表示される  
A. Java 9 以降では `-XX:MaxPermSize` が廃止されています。環境変数で `JAVA_TOOL_OPTIONS=-XX:+IgnoreUnrecognizedVMOptions` を付与してください。

---

詳細な使い方やパラメータは `config.py` や各モジュールの docstring を参照してください。  
ご質問・ご要望は Issue または開発者 (谷口: taniguchidev33@gmail.com) までご連絡ください。

