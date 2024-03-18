./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --hwe 1e-15 --keep ids/amr_ids.txt --out array_data/aou_amr_geno_hwe_maf \
  --indep-pairwise 100kb 1 0.1 --memory 100000
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep ids/amr_ids.txt \
  --exclude array_data/aou_amr_geno_hwe_maf.prune.out \
  --make-pgen --out array_data/aou_amr_geno_hwe_maf_array_pruned
./plink2 --pfile array_data/aou_amr_geno_hwe_maf_array_pruned --pca 20 approx \
  --out array_data/aou_amr_geno_hwe_maf_array_pruned
