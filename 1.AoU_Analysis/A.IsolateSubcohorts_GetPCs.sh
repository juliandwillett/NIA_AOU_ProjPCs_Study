# Use R code to process survey data to produce txt files with the IDs of individuals identifying as HISP and NONHISP, using genetics data to get individuals who identify as AMR (already determined in other study).

# First, back up the IDs into the repository
bucket="gs://fc-secure-33ab7182-9bff-4fee-bf6d-9e3b5f072033" # the Hisp/AMR study workspace
gsutil -m cp -rn ids/* $bucket/data/ids/ # back up id files (for IDs of individuals who identify as Hispanic/nonHispanic or have AMR or EUR genetic ancestry)

# Second, get the PCs for subcohorts, then back up this data. Will operate on array data, given that analysis will be similarly run in Regenie
gsutil -u $GOOGLE_PROJECT -m cp -r gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/array_data_for_regenie_step1/arrays_autosomes_post_qc_pruned_common* array_data/

# Get the PCs by the group studied. Array data already filtered for individuals who also gave srWGS
#grps=(amr eur hisp_identifying non_hisp_identifying) # not doing amr and eur, because already done previously
grps=(hisp_identifying non_hisp_identifying) ;\
for grp in "${grps[@]}"; do \
  ../plink2 --pfile array_data/arrays_autosomes_post_qc_pruned_common --pca 20 approx \
    --keep ids/${grp}_ids.txt --out pcs/${grps}_pcs ;\
done

# Backup results
gsutil -m cp -rn pcs/*eigenval $bucket/data/pc_data/
