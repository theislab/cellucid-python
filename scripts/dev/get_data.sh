#!/usr/bin/env bash
set -euo pipefail

############################################################
# USER CONFIG (LOCAL)
############################################################
# Path to file containing password for sshpass
password_file="${password_file:-/path/to/password.txt}"

# Where to put everything on your LOCAL machine
local_target_dir="/Users/kemalinecik/git_nosync/_/cellucid-data/data/raw"

# Remote login
remote_user="kemal.inecik"
remote_host="hpc-build01"

############################################################
# REMOTE DIRECTORIES
############################################################
remote_sctram_data_dir="/lustre/groups/ml01/workspace/kemal.inecik/sctram_data"
remote_processed_dir="/lustre/groups/ml01/workspace/kemal.inecik/tardis_data/processed"

############################################################
# REMOTE FINAL AnnData FILES (.h5ad)
############################################################
remote_final_files=(
  "$remote_sctram_data_dir/suo_developmental_complete.h5ad"
  "$remote_sctram_data_dir/braun_developmental_complete.h5ad"
  "$remote_sctram_data_dir/garcia_developmental_complete.h5ad"
  "$remote_sctram_data_dir/he_developmental_complete.h5ad"
  "$remote_sctram_data_dir/kanemaru_developmental_complete.h5ad"
  "$remote_sctram_data_dir/miller_developmental_complete.h5ad"
  "$remote_sctram_data_dir/norman_sciplex_cpa.h5ad"
)

############################################################
# REMOTE VARNAMES PICKLES FOR RELEVANT DATASETS
# (from your listing)
############################################################
remote_pickle_files=(
  "$remote_processed_dir/dataset_complete_Suo_varnames.pickle"
  "$remote_processed_dir/dataset_complete_Braun_varnames.pickle"
  "$remote_processed_dir/dataset_complete_Garcia_varnames.pickle"
  "$remote_processed_dir/dataset_complete_He_varnames.pickle"
  "$remote_processed_dir/dataset_complete_Kanemaru_varnames.pickle"
  "$remote_processed_dir/dataset_complete_Miller_varnames.pickle"
)

############################################################
# ENSURE LOCAL TARGET EXISTS
############################################################
mkdir -p "$local_target_dir"

############################################################
# STEP 1 — SHOW REMOTE SIZES (AnnData + varnames)
############################################################
echo "==========================================="
echo " Checking remote file sizes on $remote_host"
echo "==========================================="

remote_cmd=$(
cat <<EOF
echo "===== REMOTE FINAL AnnData SIZES (.h5ad) ====="
for f in ${remote_final_files[*]}; do
  if [ -f "\$f" ]; then
    du -h "\$f"
  else
    echo "MISSING: \$f"
  fi
done

echo ""
echo "===== REMOTE VARNAMES PICKLES SIZES (.pickle) ====="
for f in ${remote_pickle_files[*]}; do
  if [ -f "\$f" ]; then
    du -h "\$f"
  else
    echo "MISSING: \$f"
  fi
done

echo ""
echo "===== REMOTE TOTAL SIZE (AnnData + varnames) ====="
du -ch ${remote_final_files[*]} ${remote_pickle_files[*]} 2>/dev/null | grep total || true
echo "==========================================="
EOF
)

sshpass -f "$password_file" ssh -t -o LogLevel=error \
  "$remote_user@$remote_host" "$remote_cmd"

############################################################
# STEP 2 — COPY FINAL .h5ad FILES (SKIP IF ALREADY LOCAL)
############################################################
echo ""
echo "==========================================="
echo " Copying FINAL .h5ad files to local (skip existing):"
echo "   $local_target_dir"
echo "==========================================="

for f in "${remote_final_files[@]}"; do
  base_name=$(basename "$f")
  local_path="$local_target_dir/$base_name"

  echo ""
  echo "-------------------------------------------"
  echo " AnnData: $f"
  echo " Local:   $local_path"
  echo "-------------------------------------------"

  if [ -f "$local_path" ]; then
    echo "SKIP: local file already exists."
    continue
  fi

  sshpass -f "$password_file" rsync -avh --progress --ignore-existing \
    "$remote_user@$remote_host:$f" \
    "$local_target_dir/" || echo "WARN: rsync failed for $f"
done

############################################################
# STEP 3 — COPY VARNAMES PICKLES (SKIP IF ALREADY LOCAL)
############################################################
echo ""
echo "==========================================="
echo " Copying varnames pickles to local (skip existing):"
echo "   $local_target_dir"
echo "==========================================="

for f in "${remote_pickle_files[@]}"; do
  base_name=$(basename "$f")
  local_path="$local_target_dir/$base_name"

  echo ""
  echo "-------------------------------------------"
  echo " Varnames: $f"
  echo " Local:    $local_path"
  echo "-------------------------------------------"

  if [ -f "$local_path" ]; then
    echo "SKIP: local file already exists."
    continue
  fi

  sshpass -f "$password_file" rsync -avh --progress --ignore-existing \
    "$remote_user@$remote_host:$f" \
    "$local_target_dir/" || echo "WARN: rsync failed for $f"
done

############################################################
# STEP 4 — SHOW LOCAL SIZES AFTER SYNC
############################################################
echo ""
echo "==========================================="
echo " Local sizes in: $local_target_dir"
echo "==========================================="

all_remote_files=( "${remote_final_files[@]}" "${remote_pickle_files[@]}" )

for f in "${all_remote_files[@]}"; do
  base_name=$(basename "$f")
  local_path="$local_target_dir/$base_name"
  if [ -f "$local_path" ]; then
    du -h "$local_path"
  else
    echo "LOCAL MISSING: $local_path"
  fi
done

echo ""
echo "----------- LOCAL TOTAL SIZE -------------"
local_existing=()
for f in "${all_remote_files[@]}"; do
  base_name=$(basename "$f")
  lp="$local_target_dir/$base_name"
  if [ -f "$lp" ]; then
    local_existing+=( "$lp" )
  fi
done

if [ "${#local_existing[@]}" -gt 0 ]; then
  du -ch "${local_existing[@]}" 2>/dev/null | grep total || true
else
  echo "No local files found."
fi

echo "==========================================="
echo " Done."
echo "==========================================="
