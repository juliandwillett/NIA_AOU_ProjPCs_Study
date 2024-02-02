# Use R code to process survey data to produce txt files with the IDs of individuals identifying as HISP and NONHISP, using genetics data to get individuals who identify as AMR (already determined in other study).

# First, back up the IDs into the repository
bucket="gs://fc-secure-33ab7182-9bff-4fee-bf6d-9e3b5f072033" # the Hisp/AMR study workspace
gsutil -m cp -rn * $bucket/data/subcohort_ids/ # back up id files (for IDs of individuals who identify as Hispanic/nonHispanic or have AMR or EUR genetic ancestry)

# Second, get the PCs for subcohorts, then back up this data
