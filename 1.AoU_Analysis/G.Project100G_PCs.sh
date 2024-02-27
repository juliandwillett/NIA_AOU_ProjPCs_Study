# extra 1000g variants
awk 'NR==1 {print $3} NR>1 {$3=$1 "-" $2 "-" $4 "-" $5; print $3}' piezo2_work/1000G_Hg38/1000g_hg38.pvar >\
  piezo2_work/1000G_Hg38/1000g_hg38_variants.txt

./plink2 --pfile pgen_geno_1e-1_mac_20/chr18 --keep piezo2_work/rg_input/hisp_ids_3sd.txt \
   --make-bed \
  --out piezo2_work/1000G_Hg38/AoU_1000g_variants_piezo2 --chr 18 --from-bp 10644647 --to-bp 11644647

./plink --bfile piezo2_work/1000G_Hg38/AoU_1000g_variants_piezo2 --extract piezo2_work/1000G_Hg38/1000g_hg38_variants.txt \
  --make-bed --out piezo2_work/1000G_Hg38/AoU_1000g_variants_piezo2_plink19

# Do KING QC
./king -b piezo2_work/1000G_Hg38/AoU_1000g_variants_piezo2_plink19.bed --related
