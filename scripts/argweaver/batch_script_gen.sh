#!/bin/bash

# -------- user settings (edit if needed) --------
IN_TSV="scaffold_windows_150Mb.tsv"   # input with columns: scaffold  start  end
ROOT="/scratch/arghavan/LP/argweaver/script_generation"
SCRIPTS_DIR="${ROOT}/scripts"
VCF="/scratch/arghavan/LP/subsampling/s100_All_1a1b_renamed.filtered.vcf.gz"

# Slurm resources
SLURM_PART="bigmem"
SLURM_TIME="48:00:00"
SLURM_NODES=1
SLURM_MEM_PER_CPU="256G"
SLURM_NTASKS=1
SLURM_MAIL="arghavana85@gmail.com"
SLURM_MAILTYPE="FAIL"

# ARGweaver params
ARG_N=10000
ARG_R="1.6e-8"
ARG_M="1.8e-8"
ARG_NTIMES=20
ARG_MAXTIME="2e5"
ARG_C=20
ARG_NCHAINS=50   

# -------- helpers --------
sanitize() {  # keep only [A-Za-z0-9], map others to '_', squeeze repeats
  echo -n "$1" | tr -c 'A-Za-z0-9' '_' | sed -E 's/_+/_/g; s/^_|_$//g'
}
pad9() { printf "%09d" "$1"; }
randseed() {
  if command -v shuf >/dev/null 2>&1; then
    shuf -i 1-2147483647 -n 1
  else
    od -vAn -N4 -tu4 < /dev/urandom | tr -d ' '
  fi
}

mkdir -p "${ROOT}" "${SCRIPTS_DIR}"

# -------- generate one job per row --------
# Skip header if present; accept tab or space as separator.
tail -n +2 "${IN_TSV}" | awk 'NF>=3' | while read -r scaffold start end _; do
  # basic validation
  [[ -z "${scaffold:-}" || -z "${start:-}" || -z "${end:-}" ]] && continue

  safe_scaf="$(sanitize "${scaffold}")"
  tag="${safe_scaf}_$(pad9 "${start}")_$(pad9 "${end}")"
  OUT_DIR="/scratch/arghavan/LP/argweaver/final/${tag}"
  OUT_PREFIX="${OUT_DIR}/outargs_${tag}"
  LOG_OUT="${OUT_DIR}/out.txt"
  LOG_ERR="${OUT_DIR}/err.txt"
  REGION="${scaffold}:${start}-${end}"
  JOB_NAME="arg_${tag}"
  SEED=1753216431

  # mkdir -p "${OUT_DIR}"

  SH_PATH="${SCRIPTS_DIR}/${JOB_NAME}.sh"
  cat > "${SH_PATH}" <<EOF
#!/usr/bin/env bash
#SBATCH -p ${SLURM_PART}
#SBATCH -t ${SLURM_TIME}
#SBATCH -N ${SLURM_NODES}
#SBATCH --mem-per-cpu=${SLURM_MEM_PER_CPU}
#SBATCH --ntasks-per-node=${SLURM_NTASKS}
#SBATCH --mail-user=${SLURM_MAIL}
#SBATCH --mail-type=${SLURM_MAILTYPE}
#SBATCH --output="${LOG_OUT}"
#SBATCH --error="${LOG_ERR}"
#SBATCH --job-name="${JOB_NAME}"


echo "Run started at \$(date +'%Y-%m-%d %H:%M:%S')"
mkdir -p "${OUT_DIR}"

# PATHs
export PATH=\$PATH:/scratch/arghavan/mybin/argweaver/bin
export PATH=\$PATH:/scratch/arghavan/mybin/ARGweaver/bin
export PATH=\$PATH:/scratch/arghavan/mybin/ARGweaver/argweaver
export PATH=\$PATH:/scratch/arghavan/mybin/ARGweaver
export PYTHONPATH="/scratch/arghavan/mybin/ARGweaver:\$PYTHONPATH"

# 1) ARGweaver MCMC on the region
arg-sample \\
  --vcf "${VCF}" \\
  --region "${REGION}" \\
  -N ${ARG_N} --randseed 1753216431 -r ${ARG_R} -m ${ARG_M} \\
  --ntimes ${ARG_NTIMES} --maxtime ${ARG_MAXTIME} -c ${ARG_C} -n ${ARG_NCHAINS} \\
  -o "${OUT_PREFIX}"

echo "Run1 ended at \$(date +'%Y-%m-%d %H:%M:%S')"

# 2) Posterior-mean TMRCA at every site
arg-extract-tmrca "${OUT_PREFIX}.%d.smc.gz" > "${OUT_PREFIX}.tmrca.txt"

echo "Run2 at \$(date +'%Y-%m-%d %H:%M:%S')"
EOF

  chmod 0755 "${SH_PATH}"
  echo "Wrote ${SH_PATH}"
done

# -------- master submitter --------
MASTER="${SCRIPTS_DIR}/submit_all.sh"
cat > "${MASTER}" <<'EOS'
#!/usr/bin/env bash

for f in "$(dirname "$0")"/arg_*.sh; do
  echo "Submitting ${f}"
  sbatch "${f}"
done
EOS
chmod 0755 "${MASTER}"

echo "All set."
echo "Job scripts: ${SCRIPTS_DIR}"
echo "Submit them with: ${MASTER}"
