./plink2 --pfile pgen_geno_1e-1_mac_20/chr18 --keep piezo2_work/ids/hisp_identifying_ids.txt \
  --make-pgen --out piezo2_work/tmp/chr18_hisp

./plink2 --pfile piezo2_work/tmp/chr18_hisp --ld 18-11144647-G-A 18-11142027-A-G 'hwe-midp' --memory 10000

# The lead variant (rs4573993) from NIAGADS is detected in that data, why not in the GWAS data?
# Is the p value just not significant at all?

