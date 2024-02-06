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
sample_flags = list()
for (stdev in 1:5) {
    print(paste("Curr stdev:",stdev))
    sample_flags[[stdev]] = lapply(eigenvec[[1]],function(x) {
        row = eigenvec %>% filter(`#IID` == x)
        curr_flag = character()
        for (col in 2:21) {
            calc = (row[[col]] - mean_sd[[col-1]][[1]]) / (stdev * mean_sd[[col-1]][[2]])
            if (abs(calc) > 5)
                curr_flag %<>% append(glue("PC{col-1}"))
        }
        if (length(curr_flag) > 0) paste0(curr_flag,collapse = ",")
        else NA
    })
}
sample_flags_hisp = data.frame(IID=eigenvec[[1]],One_SD_Flags = t(as.data.frame(sample_flags[[1]])),
                              Two_SD_Flags = t(as.data.frame(sample_flags[[2]])),
                              Three_SD_Flags = t(as.data.frame(sample_flags[[3]])),
                              Four_SD_Flags = t(as.data.frame(sample_flags[[4]])),
                               Five_SD_Flags = t(as.data.frame(sample_flags[[5]])))
vroom_write(sample_flags_hisp,"piezo2_work/pop_qc/hisp_sample_flags_by_sd.txt")
print("wrote to file")

### Check PCs for Non HISP
eigenvec_nonhisp = vroom("piezo2_work/pcs/non_hisp_identifying_pcs.eigenvec",show_col_types=F)
mean_sd = lapply(eigenvec[2:21],FUN = function(x) {
    c(mean(x),sd(x))
})

sample_flags = list()
for (stdev in 1:5) {
    print(paste("Curr stdev:",stdev))
    sample_flags[[stdev]] = lapply(eigenvec[[1]],function(x) {
        row = eigenvec %>% filter(`#IID` == x)
        curr_flag = character()
        for (col in 2:21) {
            calc = (row[[col]] - mean_sd[[col-1]][[1]]) / (stdev * mean_sd[[col-1]][[2]])
            if (abs(calc) > 5)
                curr_flag %<>% append(glue("PC{col-1}"))
        }
        if (length(curr_flag) > 0) paste0(curr_flag,collapse = ",")
        else NA
    })
}

sample_flags_non_hisp = data.frame(IID=eigenvec[[1]],One_SD_Flags = t(as.data.frame(sample_flags[[1]])),
                              Two_SD_Flags = t(as.data.frame(sample_flags[[2]])),
                              Three_SD_Flags = t(as.data.frame(sample_flags[[3]])),
                              Four_SD_Flags = t(as.data.frame(sample_flags[[4]])),
                               Five_SD_Flags = t(as.data.frame(sample_flags[[5]])))
vroom_write(sample_flags_non_hisp,"piezo2_work/pop_qc/non_hisp_sample_flags_by_sd.txt")
print("wrote to file")

###############
# Plot the results, to see the ideal number of SD to use:
anc = vroom("ancestry_preds.tsv",show_col_types = F) %>% rename(IID = research_id)
flagged_samples = vroom("piezo2_work/pop_qc/hisp_sample_flags_by_sd.txt",show_col_types=F)
eigenvec_hisp = vroom("piezo2_work/pcs/hisp_identifying_pcs.eigenvec",show_col_types=F) %>% 
    rename(IID = `#IID`) %>% left_join(anc,by = "IID") %>% left_join(flagged_samples,by="IID")
head(flagged_samples)

pdf(file = "hisp_pc_sd_plots.pdf")
for (n in 1:5) {
    print(ggplot(eigenvec_hisp %>% filter(is.na(.[[25+n]])),aes(x=PC1,y=PC2,color=ancestry_pred)) + geom_point() +
        theme(text = element_text(size=14)) + ggtitle(glue("SD Cutoff: {n}")))
}
dev.off()

###
# So go for 3 SDs
out_df = flagged_samples %>% filter(is.na(Three_SD_Flags)) %>% select(IID) %>%
    mutate(FID = 0,.before = "IID")
vroom_write(out_df,"piezo2_work/rg_input/hisp_ids_3sd.txt")
