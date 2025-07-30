FROM continuumio/miniconda3:23.3.1-0

# 環境名
ENV ENV_NAME=auto-app2-env

# 作業ディレクトリ
WORKDIR /app

# 必要ファイルコピー
COPY environment.yml .
COPY . /app

# Conda環境作成（Docker build時に実行）
RUN conda env create -f environment.yml

# conda activate を有効にする設定
SHELL ["conda", "run", "-n", "auto-app2-env", "/bin/bash", "-c"]

# アプリケーションファイルをコピー
COPY . /app

# データディレクトリの作成
RUN mkdir -p /app/data/logs /app/data/raw_data /app/data/reference /app/data/results /app/data/temp

# 権限の設定
RUN chmod +x /app/main.py

# 環境変数の設定
ENV PATH="/opt/conda/envs/auto-app2-env/bin:$PATH"

# デフォルトコマンド（bashで起動）
CMD ["conda", "run", "-n", "auto-app2-env", "bash"]