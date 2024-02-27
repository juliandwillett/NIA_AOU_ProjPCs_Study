# Get subsetted plink data
./plink2 --pfile /n/holystore01/LABS/tanzi_lab/Users/dmitry/NIAGADS_v9_analysis/NIAGADS_pgen/cc34438_combined \
  --keep /n/holystore01/LABS/tanzi_lab/Users/dmitry/NIAGADS_v9_analysis/keep_8467_selfHISP.txt \
  --chr 18 --from-bp 10644647 --to-bp 11644647 --make-pgen --out \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_PIEZO2_PLINK --memory 80000

# Get LD information. Variants in format of: 18:10723405:G:A
