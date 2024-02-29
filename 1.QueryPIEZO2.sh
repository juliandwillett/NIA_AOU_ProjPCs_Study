# Get subsetted plink data
./plink2 --pfile /n/holystore01/LABS/tanzi_lab/Users/dmitry/NIAGADS_v9_analysis/NIAGADS_pgen/cc34438_combined \
  --keep /n/holystore01/LABS/tanzi_lab/Users/dmitry/NIAGADS_v9_analysis/keep_8467_selfHISP.txt \
  --chr 18 --from-bp 10644647 --to-bp 11644647 --make-pgen --out \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_PIEZO2_PLINK --memory 80000

# Get LD information. Variants in format of: 18:10723405:G:A
./plink2 --pfile /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_PIEZO2_PLINK \
  --r2-phased --ld-snp-list /n/home09/jwillett/true_lab_storage/02_PIEZO2_HISP/for_niagads_ld.txt \
  --out /n/home09/jwillett/true_lab_storage/02_PIEZO2_HISP/piezo2_ld --ld-window-r2 0

#######
# Get info for fine mapping: 500 kb window
# Get all variants in window, no p value filter
zcat /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.dn8.gz | \
awk 'NR==1 || ($1 == 18 && $2 >= 10894647 && $2 < 11394647 && $8 != "NA" && $9 != "NA") {print}' > \
/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.piezo2.txt

# Get the IDs
awk 'NR>1 && $8 != "NA" && $9 != "NA" {print $3}' \
/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.piezo2.txt > \
/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.piezo2_just_ids.txt

# Get the R matrix for susieR
./plink2 --pfile /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_PIEZO2_PLINK \
  --r-phased square --extract /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.piezo2_just_ids.txt \
  --out /n/home09/jwillett/true_lab_storage/02_PIEZO2_HISP/piezo2_finemap_data
