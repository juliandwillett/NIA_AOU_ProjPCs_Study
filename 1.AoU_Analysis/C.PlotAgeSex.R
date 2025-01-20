# Get data. Note that Female Sex in NIAGADS was 1, and 0 in AoU. So swap as necessary.
pheno_df = vroom('regenie_input/regenie_pheno.txt',show_col_types = F)
nia_hisp_cohort = vroom("NIAGADS_demographic_data_all.txt",show_col_types = F) %>%
  mutate(Sex = ifelse(Sex == 1,0,1))
aou_hisp_cohort = merge(vroom("regenie_input/regenie_covar_hisp_anc_covar.txt",show_col_types = F),
                       pheno_df,by='IID')
aou_amr_cohort = merge(vroom("regenie_input/regenie_covar_commonpcs.txt",show_col_types = F) %>%
                          filter(IID %in% vroom("amr_ids.txt",show_col_types = F,col_names = F)$X2),
                      pheno_df,by='IID')
aou_nia_proj_3sd = merge(vroom("regenie_input/regenie_covar_projected_niagads_hisp_3sd.txt",show_col_types = F),
                        pheno_df,by='IID')
aou_nia_proj_matchit = merge(vroom("regenie_input/aou_nia_matchit_noanccovar.txt",show_col_types = F),
                            pheno_df,by='IID')

# Put in a single data frame for plotting
all_cohort_df = data.frame(Age=nia_hisp_cohort$Age,Sex=nia_hisp_cohort$Sex,AD=nia_hisp_cohort$Affection.Status,Cohort='NIAGADS HISP') %>%
    add_row(data.frame(Age = aou_hisp_cohort$Age,Sex = aou_hisp_cohort$Sex,AD=aou_hisp_cohort$AD_any,Cohort = 'AoU HISP')) %>%
    add_row(data.frame(Age = aou_amr_cohort$Age,Sex = aou_amr_cohort$Sex,AD=aou_amr_cohort$AD_any,Cohort = 'AoU AMR')) %>%
    add_row(data.frame(Age = aou_nia_proj_3sd$Age,Sex = aou_nia_proj_3sd$Sex,AD=aou_nia_proj_3sd$AD_any,Cohort = 'AoU NIA Proj 3SD')) %>%
    add_row(data.frame(Age = aou_nia_proj_matchit$Age,Sex = aou_nia_proj_matchit$Sex,AD=aou_nia_proj_matchit$AD_any,Cohort = 'AoU NIA Proj MatchIt'))

# Plot and save
# AGE PLOT
cud_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plt = ggplot(all_cohort_df %>% mutate(AD = as.factor(AD),AD = ifelse(AD == "1","Case","Control")),
       aes(x=Cohort,y=Age,fill=AD)) + geom_boxplot(position=position_dodge(1)) +
    theme_bw() + theme(text = element_text(size=20)) + scale_fill_manual(values = cud_palette)
print(plt)
ggsave(plt,'age_plot.png',width = 5,height=3)

# SEX PLOT. 
cud_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
proportions <- all_cohort_df %>% group_by(Cohort, AD) %>%
  summarise(Proportion = mean(Sex == 0))
plt = ggplot(proportions %>% mutate(AD = as.factor(AD),AD = ifelse(AD == "1","Case","Control")), 
             aes(x = Cohort, y = Proportion, fill = AD)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(round(Proportion * 100, 0), "%")), position = position_dodge(width = 0.9), vjust = -0.5) +
  theme_bw() + theme(text = element_text(size=20)) + scale_fill_manual(values = cud_palette) +
    ylim(c(0,0.8))                  
print(plt)
ggsave("sex_plot.png", plot = plt, width = 12, height = 4, dpi = 300)
