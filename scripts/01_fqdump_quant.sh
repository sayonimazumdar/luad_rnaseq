#!/usr/bin/env bash
# 01_fqdump_quant.sh : per-SAMPLE loop using NCBI fastq-dump -X (reads only the first N reads of each run over the network).
# Steps per sample: fastq-dump each run (with retries) -> fastp -> ONE Salmon quant -> delete reads.
# Run from the project root with the luad-rnaseq env active:  bash scripts/01_fqdump_quant.sh
set -euo pipefail

THREADS=4
MAX_SPOTS=6000000                  # TOTAL reads per sample, split evenly across its runs
INDEX="ref/salmon_index"
SAMPLES="metadata/samples.csv"

mkdir -p fastq trimmed quants qc logs

# download the first N read(s) of a run; retry if the connection drops or the mates don't match
dump_run () {   # run n layout
    local run="$1" n="$2" layout="$3" try l1 l2
    for try in 1 2 3 4 5; do
        rm -f fastq/${run}*
        if fastq-dump -X "${n}" --split-files --skip-technical \
               --defline-seq '@$ac.$si/$ri' --defline-qual '+' -O fastq "${run}" 2> "logs/${run}_fqdump.log"; then
            if [[ "${layout}" == "PAIRED" ]]; then
                if [[ -s "fastq/${run}_1.fastq" && -s "fastq/${run}_2.fastq" ]]; then
                    l1=$(wc -l < "fastq/${run}_1.fastq"); l2=$(wc -l < "fastq/${run}_2.fastq")
                    [[ "${l1}" -eq "${l2}" ]] && return 0
                fi
            else
                [[ -s "fastq/${run}.fastq" ]] && return 0
            fi
        fi
        echo "  fastq-dump for ${run} failed (try ${try}/5), retrying in 30s" >&2
        sleep 30
    done
    return 1
}

SAMPLE_IDS=$(tail -n +2 "${SAMPLES}" | tr -d '\r' | awk -F, '{print $2"_"$3}' | sort -u)

for SID in ${SAMPLE_IDS}; do
    if [[ -s "quants/${SID}/quant.sf" ]]; then
        echo "[$(date +%T)] ${SID}: already quantified, skipping"; continue
    fi

    RUNS=$(tail -n +2 "${SAMPLES}" | tr -d '\r' | awk -F, -v s="${SID}" '$2"_"$3==s {print $1}')
    LAYOUT=$(tail -n +2 "${SAMPLES}" | tr -d '\r' | awk -F, -v s="${SID}" '$2"_"$3==s {print toupper($4); exit}')
    NRUNS=$(echo "${RUNS}" | wc -w)
    PER_RUN=$(( MAX_SPOTS / NRUNS ))
    echo "[$(date +%T)] ==== ${SID}: ${NRUNS} run(s) [$(echo ${RUNS})], ${LAYOUT}, ${PER_RUN} reads per run ===="

    R1=(); R2=()
    for RUN in ${RUNS}; do
        dump_run "${RUN}" "${PER_RUN}" "${LAYOUT}" || { echo "fastq-dump kept failing for ${RUN}"; exit 1; }
        if [[ "${LAYOUT}" == "PAIRED" ]]; then
            fastp -i "fastq/${RUN}_1.fastq" -I "fastq/${RUN}_2.fastq" \
                  -o "trimmed/${RUN}_1.fastq.gz" -O "trimmed/${RUN}_2.fastq.gz" \
                  -w "${THREADS}" -j "qc/${RUN}_fastp.json" -h "qc/${RUN}_fastp.html" 2> "logs/${RUN}_fastp.log"
            R1+=("trimmed/${RUN}_1.fastq.gz"); R2+=("trimmed/${RUN}_2.fastq.gz")
        else
            fastp -i "fastq/${RUN}.fastq" -o "trimmed/${RUN}.fastq.gz" \
                  -w "${THREADS}" -j "qc/${RUN}_fastp.json" -h "qc/${RUN}_fastp.html" 2> "logs/${RUN}_fastp.log"
            R1+=("trimmed/${RUN}.fastq.gz")
        fi
        rm -f fastq/${RUN}*
    done

    if [[ "${LAYOUT}" == "PAIRED" ]]; then
        salmon quant -i "${INDEX}" -l A -1 "${R1[@]}" -2 "${R2[@]}" -p "${THREADS}" --validateMappings \
            -o "quants/${SID}" 2> "logs/${SID}_salmon.log"
    else
        salmon quant -i "${INDEX}" -l A -r "${R1[@]}" -p "${THREADS}" --validateMappings \
            -o "quants/${SID}" 2> "logs/${SID}_salmon.log"
    fi
    rm -f trimmed/*
    echo "[$(date +%T)] ${SID}: done. Disk free: $(df -h . | awk 'NR==2{print $4}')"
done
echo "All samples processed."
