# Requires 70 GB of ram, takes less than an hour

wc -l no_missing_covar2 # right covar file for inclusion
awk '{gsub(/\.[0]+$/, "", $1); gsub(/\.[0]+$/, "", $2); print $1 "\t" $2}' no_missing_covar2 > fid_iid_participants.txt

./plink2 --pfile merged_chromosome_18 --keep fid_iid_participants.txt \
  --geno 0.1 \
  --r2-phased --ld-window-r2 0.1 --chr 18 --from-bp 10644647 --to-bp 11644647 \
  --out ukb_all_piezo2_geno_1e-1_r2_1e-1 --memory 100000
  
