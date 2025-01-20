setwd("~/true_lab_storage/02_NIAGADS_HISP")
library(vroom)
library(tidyverse)
library(glue)
library(magrittr)
library(writexl)
library(susieR)
source("Rcode/Functions.R")
`%notin%` = negate(`%in%`)

############
# Organize GWAS data to put results from each file into one unified file
# organize_gwas_data('MatchIt_Gen')
# organize_gwas_data('MatchIt_GenPhen')

############
# Organize GW hits from AoU, to query them for HWE and sample missingness
# gw_data_aou_nia_proj = get_gwas_nominal_hits("AoU_NIA_PROJ",p_cutoff = 5e-8)
# vroom_write(data.frame(ID=gw_data_aou_nia_proj$ID),"gw_hits_aou_nia_proj.txt",col_names = F)

###########
# Get lambda values to evaluate for genomic inflation
# lambda_aou_nia_proj = getLambda("AoU_NIA_PROJ")
# lambda_niagads = getLambda("NIAGADS_HISP") 
# saveRDS(list(lambda_aou_hisp_hisp_pcs,lambda_aou_hisp_all_pcs,lambda_aou_nia_proj,
#              lambda_aou_amr,lambda_niagads),"lambda_gc_and_1000_values.rds")

##########
# FIGURE 1: STUDY DESIGN.

##########
# FIGURE 2: AGE AND SEX STATISTICS. PRODUCED IN AOU

##########
# FIGURE 3: MANHATTAN FOR EACH GWAS
# GET HWE DATA FOR QC, following:
#   https://github.com/juliandwillett/AoU_AlzheimersGWAS/blob/main/1.AoU_GWAS/K.CheckForFalsePositives_Anc.sh
# make_files_for_hwe()

# 3A: NIAGADS HISP
gw_data_nia_hisp = get_gwas_nominal_hits("NIAGADS_HISP",p_cutoff = 1e-7)
lead_variants_nia_hisp = get_lead_variants(gw_data_nia_hisp,get_gene_names = T)
figure3a = make_manhattan_easy_label(dataset = "NIAGADS_HISP",lead_var=lead_variants_nia_hisp,
                                     fig_num="3A",p_data_cut=1e-7,p_label_cut=1e-7,show_highlights=F)

# 3B: AOU MATCHED NIA HISP GEN PCS ONLY
gw_data_aou_genpcs = get_gwas_nominal_hits("AoU_NIA_Gen_GW",p_cutoff = 5e-8)
lead_variants_aou_genpcmatch = get_lead_variants(gw_data_aou_genpcs,get_gene_names = T)
figure3b = make_manhattan_easy_label(dataset = "AoU_NIA_Gen_GW",gw_data_aou_genpcs,
                                    lead_variants_aou_genpcmatch,fig_num = "3B",show_gene_names = T,
                                    only_new_old="none",
                                    p_data_cut = 5e-8,p_label_cut = 5e-8,show_highlights = F)

# 3C: AOU MATCHED NIA HISP GEN PCS AGE SEX
gw_data_aou_genpcs_phen = get_gwas_nominal_hits("AoU_NIA_GenPhen_GW",p_cutoff = 5e-8)
lead_variants_aou_genphenmatch = get_lead_variants(gw_data_aou_genpcs_phen,get_gene_names = T)
figure3c = make_manhattan_easy_label(dataset = "AoU_NIA_GenPhen_GW",gw_data_aou_genpcs_phen,
                                     lead_variants_aou_genphenmatch,fig_num = "3C",show_gene_names = T,
                                     only_new_old="none",
                                     p_data_cut = 5e-8,p_label_cut = 5e-8,show_highlights = F)
############
# TABLE 1: LEAD VARIANTS FOR EACH STUDY
# input data
gw_data_nia_hisp = get_gwas_nominal_hits("NIAGADS_HISP",p_cutoff = 1e-7)
lead_var_nia_hisp = get_lead_variants(gw_data_nia_hisp,get_gene_names = T) %>% 
  rename(Rsid=Rsid.x,ALLELE1 = effectallele)
gw_data_aou_genpcs = get_gwas_nominal_hits("AoU_NIA_Gen_GW",p_cutoff = 5e-8)
lead_var_aou_genpcmatch = get_lead_variants(gw_data_aou_genpcs,get_gene_names = T)
gw_data_aou_genpcs_phen = get_gwas_nominal_hits("AoU_NIA_GenPhen_GW",p_cutoff = 5e-8)
lead_var_aou_genphenmatch = get_lead_variants(gw_data_aou_genpcs_phen,get_gene_names = T)
# intersect_for_p_values(lead_var_nia_hisp,lead_var_aou_genpcmatch,lead_var_aou_genphenmatch)
# table function
table1 = make_gwas_lead_locus_table(lead_var_nia_hisp,lead_var_aou_genpcmatch,lead_var_aou_genphenmatch)
write_xlsx(table1,"Paper_Tables_Figures/table1.xlsx")

############
# FIGURE 4: MANHATTANS FOR META ANALYSES
# Fig 4a - GEN PC only match
meta_data_genpc_match = get_meta_nominal_hits("AoU_NIA_Gen_GW",5e-8) %>% filter(HetDf > 0)
lead_var_meta_genpc_match = get_lead_variants(meta_data_genpc_match,get_gene_names = T,meta = T)
figure4a = make_manhattan_easy_label(dataset = "Meta_AoU_NIA_Gen_GW",meta_data_genpc_match,
                                     lead_var_meta_genpc_match,fig_num = "4A",show_gene_names = T,
                                     only_new_old="none",
                                     p_data_cut = 5e-8,p_label_cut = 5e-8,show_highlights = F)

# Fig 4b - GEN PC and Age and Sex match
meta_data_genpc_phen_match = get_meta_nominal_hits("AoU_NIA_GenPhen_GW",5e-8) %>% filter(HetDf > 0)
lead_var_meta_genpc_phen_match = get_lead_variants(meta_data_genpc_phen_match,get_gene_names = T,meta=T)
figure4b = make_manhattan_easy_label(dataset = "Meta_AoU_NIA_GenPhen_GW",meta_data_genpc_phen_match,
                                     lead_var_meta_genpc_phen_match,fig_num = "4B",show_gene_names = T,
                                     only_new_old="none",
                                     p_data_cut = 5e-8,p_label_cut = 5e-8,show_highlights = F)

##################
# TABLE 2 LEAD VARIANT INFO FOR META
# get necessary statistics
meta_data_genpc_match = get_meta_nominal_hits("AoU_NIA_Gen_GW",5e-8) %>% filter(HetDf > 0)
lead_var_meta_genpc_match = get_lead_variants(meta_data_genpc_match,get_gene_names = T,meta = T)
meta_data_genpc_phen_match = get_meta_nominal_hits("AoU_NIA_GenPhen_GW",5e-8) %>% filter(HetDf > 0)
lead_var_meta_genpc_phen_match = get_lead_variants(meta_data_genpc_phen_match,get_gene_names = T,meta=T)
# make file for intersection (for p vals)
# vroom_write(data.frame(CHRPOS=meta_data_genpc_match$CHRPOS) %>%
#               add_row(data.frame(CHRPOS=meta_data_genpc_phen_match$CHRPOS)),
#             "working/chrpos_for_meta_intersect.txt")
# system("bash intersect_meta.sh")
# make table
table2 = make_meta_lead_locus_table(lead_var_meta_genpc_match,lead_var_meta_genpc_phen_match)
write_xlsx(table2,"Paper_Tables_Figures/table2.xlsx")

#################
# TABLE 3: GENE GENERAL ANNOTATION:
data = table2 %>% filter(NewOld == "New",NIA_P<=0.05,AoU_Gen_P < 0.05 | AoU_GenPhen_P < 0.05)
table3 = get_gene_expression(data)
write_xlsx(table3,"Paper_Tables_Figures/table3.xlsx")

#################
# TABLE 4: GENE PATH SPECIFIC ANNOTATION
table4 = get_pathspecific_gene_expression(data)
write_xlsx(table4,"Paper_Tables_Figures/table4.xlsx")
