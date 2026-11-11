#!/usr/bin/env bash
# Background fetches for assets that are needed before E1' can run.
#   - bpRNA-1m(90): the large training pool (structure labels included)
#   - RiNALMo-giga weights (AGPL-3.0): frozen backbone
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
mkdir -p "${ROOT}/data/raw/bprna" "${ROOT}/data/models"

if [ ! -s "${ROOT}/data/raw/bprna/bpRNA_1m.json" ]; then
  echo "[bprna] downloading bpRNA-1m from HuggingFace mirror"
  for base in "https://huggingface.co" "https://hf-mirror.com"; do
    wget -c --tries=5 --waitretry=10 --timeout=120 \
      -O "${ROOT}/data/raw/bprna/bpRNA_1m.json" \
      "${base}/datasets/rouskinlab/bpRNA-1m/resolve/main/data.json" && break
  done
fi

if [ ! -x "${ROOT}/tools/mmseqs/bin/mmseqs" ]; then
  echo "[mmseqs] downloading static build"
  mkdir -p "${ROOT}/tools/mmseqs"
  wget -c --tries=10 --waitretry=10 --timeout=120 \
    -O "${ROOT}/tools/mmseqs.tar.gz" \
    "https://github.com/soedinglab/MMseqs2/releases/download/15-6f452/mmseqs-linux-avx2.tar.gz" \
    && tar xzf "${ROOT}/tools/mmseqs.tar.gz" -C "${ROOT}/tools" \
    && rm -rf "${ROOT}/tools/mmseqs" \
    && mv "${ROOT}/tools/mmseqs-linux-avx2" "${ROOT}/tools/mmseqs"
fi

if [ ! -d "${ROOT}/data/models/rinalmo-giga" ]; then
  echo "[rinalmo] downloading giga"
  HF_HOME="${ROOT}/data/models/.hf" "${PY}" - <<'PY'
from huggingface_hub import snapshot_download
p = snapshot_download(
    "multimolecule/rinalmo-giga",
    local_dir="data/models/rinalmo-giga",
    allow_patterns=["*.json", "*.safetensors", "*.txt", "*.md", "*.py"],
)
print("downloaded to", p)
PY
fi

echo "[done]"
ls -la "${ROOT}/data/raw/bprna" "${ROOT}/data/models" 2>/dev/null
