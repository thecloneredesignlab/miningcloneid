#!/usr/bin/env bash
set -euo pipefail

echo "[download.sh] Downloading DetectionResults_InVitro..."
mkdir -p data

if ! command -v aws >/dev/null; then
  echo "ERROR: aws CLI not found"; exit 1
fi

if [ -z "${AWS_ACCESS_KEY_ID:-}" ] || [ -z "${AWS_SECRET_ACCESS_KEY:-}" ]; then
  echo "WARNING: AWS credentials not set; skipping download."; exit 0
fi

aws s3 cp "s3://cloneid4mysql8/DetectionResults_InVitro.tar.gz" - | tar -xz -C data

echo "[download.sh] Done."
