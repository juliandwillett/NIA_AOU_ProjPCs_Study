# get files for this process:
bucket="gs://fc-secure-33ab7182-9bff-4fee-bf6d-9e3b5f072033" # the Hisp/AMR study workspace
gsutil -m cp -rn $bucket/data/pc_data/*eigenvec pcs/*eigenval 

{R}
### Check PCs for HISP
# Get the mean and SD
eigenvec = vroom("piezo2_work/pcs/hisp_identifying_pcs.eigenvec",show_col_types=F)
mean_sd = lapply(eigenvec[2:21],FUN = function(x) {
    c(mean(x),sd(x))
})

# Check for samples that go beyond 5 SDs
sample_flags = lapply(eigenvec[[1]],function(x) {
    row = eigenvec %>% filter(`#IID` == x)
    curr_flag = character()
    for (col in 2:21) {
        calc = (row[[col]] - mean_sd[[col-1]][[1]]) / (5 * mean_sd[[col-1]][[2]])
        if (abs(calc) > 5)
            curr_flag %<>% append(glue("PC{col-1}"))
    }
    if (length(curr_flag) > 0) paste0(curr_flag,collapse = ",")
    else NA
})

num_flagged_samples = length(which(!is.na(sample_flags)))
print(glue("Num samples: {length(sample_flags)}. Num flagged samples using 20 PCs: {num_flagged_samples}"))
#df = data.frame(IID = eigenvec[[1]],Flags = sample_flags)

### Check PCs for Non HISP
eigenvec = vroom("piezo2_work/pcs/non_hisp_identifying_pcs.eigenvec",show_col_types=F)
mean_sd = lapply(eigenvec[2:21],FUN = function(x) {
    c(mean(x),sd(x))
})

# Check for samples that go beyond 5 SDs
sample_flags = lapply(eigenvec[[1]],function(x) {
    row = eigenvec %>% filter(`#IID` == x)
    curr_flag = character()
    for (col in 2:21) {
        calc = (row[[col]] - mean_sd[[col-1]][[1]]) / (5 * mean_sd[[col-1]][[2]])
        if (abs(calc) > 5)
            curr_flag %<>% append(glue("PC{col-1}"))
    }
    if (length(curr_flag) > 0) paste0(curr_flag,collapse = ",")
    else NA
})

num_flagged_samples = length(which(!is.na(sample_flags)))
print(glue("Num samples: {nrow(eigenvec)}. Num flagged samples using 20 PCs: {num_flagged_samples}"))
#df = data.frame(IID = eigenvec[[1]],Flags = sample_flags)
