# get files for this process:
bucket="gs://fc-secure-33ab7182-9bff-4fee-bf6d-9e3b5f072033" # the Hisp/AMR study workspace
gsutil -m cp -rn $bucket/data/pc_data/*eigenvec pcs/*eigenval 
