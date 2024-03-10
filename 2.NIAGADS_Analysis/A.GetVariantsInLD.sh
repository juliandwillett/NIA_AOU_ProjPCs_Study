# key variants: rs5018879, rs3925017, rs4573993, rs4413029
salloc -p test --mem 160000 -t 0-01:00 -c 6
./plink2 --pfile /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_PIEZO2_PLINK \
  --r2-phased --ld-window-r2 0.9 --out /n/home09/jwillett/true_lab_storage/02_PIEZO2_HISP/piezo2_r2_phased_matrix \
  --memory 160000
