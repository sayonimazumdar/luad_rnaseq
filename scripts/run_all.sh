#!/usr/bin/env bash
# run_all.sh : full pipeline in order (run from the project root)
# NOTE: step 01 takes many hours; run it with nohup (see README).
set -euo pipefail
bash scripts/00_setup.sh            # once
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate luad-rnaseq
bash scripts/01_fqdump_quant.sh
bash scripts/02_multiqc_tx2gene.sh
