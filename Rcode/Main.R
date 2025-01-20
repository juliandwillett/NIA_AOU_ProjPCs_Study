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
# organize_gwas_data("amr")
# organize_gwas_data("Proj_3SD")
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
# TABLE 3: FAVOR ANNOTATION: solely on variants that occur in both studies.
# NOT VERY INTERESTING, BUT ALSO OF LIMITED TRANSLATABILITY GIVEN SAMPLE SIZE
# meta_data_favor_annot = get_favor_annot(meta_data_genpc_match,meta_data_genpc_phen_match)
# table3 = make_favor_functional_table(meta_data_favor_annot)

#################
# TABLE 3: GENE GENERAL ANNOTATION:
data = table2 %>% filter(NewOld == "New",NIA_P<=0.05,AoU_Gen_P < 0.05 | AoU_GenPhen_P < 0.05)
table3 = get_gene_expression(data)
write_xlsx(table3,"Paper_Tables_Figures/table3.xlsx")

#################
# TABLE 4: GENE PATH SPECIFIC ANNOTATION
table4 = get_pathspecific_gene_expression(data)
write_xlsx(table4,"Paper_Tables_Figures/table4.xlsx")







#############
# Figure 2: demographic information. Produced in AoU.
figure_2_data = write_demographic_data_for_export()

#############
# Figure 3: Produce Manhattan plot for NIAGADS Hispanic


############
# Table 1 - locus information
### First, give lead loci
### Then, give top variants in PIEZO2 locus
gw_data_nia_hisp = get_gwas_nominal_hits("NIAGADS_HISP",p_cutoff = 1e-5)
lead_variants_nia_hisp = get_lead_variants(gw_data_nia_hisp,get_gene_names=F)
table1 = make_lead_locus_table(gw_data_nia_hisp,lead_variants_nia_hisp)
write_xlsx(table_1,"Paper_Tables_Figures/table1.xlsx")

##########
# Figure 4C, i.e. Manhattan for AoU NIA GenPhenSim Matchit. GW Hits not nominal in NIA
gw_data_aou_proj_matchit = get_gwas_nominal_hits("AoU_NIA_MatchIt_GenPhen",p_cutoff = 1e-5) 
lead_variants_aou_proj_matchit = get_lead_variants(gw_data_aou_proj_matchit,get_gene_names=T)
figure4d = make_manhattan_easy_label(dataset = "AoU_NIA_MatchIt_GenPhen",gw_data_aou_proj_matchit,
                                     lead_variants_aou_proj_matchit,NA,
                                     fig_num="4D",show_gene_names=T,p_data_cut = 1e-5,
                                     p_label_cut=1e-7)
loci = getNumberCommonRareLoci(aou_hisp_nominal_nia,maf_cut=0.01,return_lead=T)

##########
# Figure 4D, i.e. Manhattan for meta between AoU Nia GenPhenSim MatchIt against NIA HISP
gw_meta_aou_proj_genphen_matchit_nia_hisp = get_gwas_nominal_hits("AoU_NIA_GenPhen_Meta",p_cutoff = 1e-5)

gw_data_aou_proj_matchit = get_gwas_nominal_hits("AoU_NIA_MatchIt",p_cutoff = 1e-5) 
lead_variants_aou_proj_matchit = get_lead_variants(gw_data_aou_proj_matchit,get_gene_names=T)
figure4d = make_manhattan_easy_label(dataset = "AoU_NIA_MatchIt",gw_data_aou_proj_matchit,
                                     lead_variants_aou_proj_matchit,NA,
                                     fig_num="4D",show_gene_names=T,p_data_cut = 1e-5,
                                     p_label_cut=1e-7)
loci = getNumberCommonRareLoci(aou_hisp_nominal_nia,maf_cut=0.01,return_lead=T)

##########
### Produce Table 2, i.e. lead variants from AoU studies, nominal in AoU for only 3 PC dataset that had lots of power
# First, get the CHRPOS of all GW-significant hits
lead_variants_aou_hisp = get_lead_variants(get_gwas_nominal_hits("AoU_HISP",p_cutoff = 5e-8),get_gene_names=T)
lead_variants_aou_amr = get_lead_variants(get_gwas_nominal_hits("AoU_AMR",p_cutoff = 5e-8) ,get_gene_names=T)
lead_variants_aou_nia_proj = get_lead_variants(get_nia_nominal_hits(get_gwas_nominal_hits("AoU_NIA_PROJ_3SD",p_cutoff = 5e-8)),
                                               get_gene_names = F) %>% rename(Rsid = Rsid.x)
lead_variants_aou_proj_matchit = get_lead_variants(get_gwas_nominal_hits("AoU_NIA_MatchIt",p_cutoff = 5e-8),get_gene_names=T)
variant_intersect = unique(lead_variants_aou_hisp$CHRPOS %>% append(lead_variants_aou_amr$CHRPOS) %>%
  append(lead_variants_aou_nia_proj$CHRPOS) %>% append(lead_variants_aou_proj_matchit$CHRPOS))
# data_intersects = get_aou_cohort_intersects(variant_intersect) ; saveRDS(data_intersects,"aou_multicohort_hits_intersect.rds")
data_intersects = readRDS("aou_multicohort_hits_intersect.rds")
table2 = make_aou_lead_variant_table(lead_variants_aou_hisp,lead_variants_aou_amr,
                                     lead_variants_aou_nia_proj,lead_variants_aou_proj_matchit,
                                     data_intersects)
write_xlsx(table2,"Paper_Tables_Figures/table2.xlsx")

########## 
# Table 3 - functional information for PIEZO2 in NIAGADS
table3 = make_functional_table_niagads(table1)
write_xlsx(table3,"Paper_Tables_Figures/table3.xlsx")

##########
# Table 4 - functional information for nominal loci in AoU
aou_hisp = get_gwas_nominal_hits("AoU_HISP",p_cutoff = 5e-8) %>% filter(NewOld=='New')
aou_amr = get_gwas_nominal_hits("AoU_AMR",p_cutoff = 5e-8) %>% filter(NewOld=='New')
aou_nia_proj = get_nia_nominal_hits(get_gwas_nominal_hits("AoU_NIA_PROJ_3SD",p_cutoff = 5e-8)) %>% filter(NewOld=='New')
aou_nia_matchit = get_gwas_nominal_hits("AoU_NIA_MatchIt",p_cutoff = 5e-8) %>% filter(NewOld=='New')
table4 = make_functional_table_aou(aou_hisp,aou_amr,aou_nia_proj,aou_nia_matchit)
write_xlsx(table4,"Paper_Tables_Figures/table4.xlsx")

#########
# Table 5 - cell expr table, all path features
sc_data = getSingleCellData(0.01)
table3 = read_xlsx("Paper_Tables_Figures/table3.xlsx")
table4 = read_xlsx("Paper_Tables_Figures/table4.xlsx")
table5 = make_general_expr_table(table3,table4,sc_data$Gross_Phenotypes)
write_xlsx(table5,"Paper_Tables_Figures/table5.xlsx")

#########
# Table 6 - cell expr table, by path feature
table6 = make_path_specific_expr_table(table3,table4,sc_data$Micro_Phenotypes)
write_xlsx(table6,"Paper_Tables_Figures/table6.xlsx")

########
# Supplemental Figure 1, i.e. Manhattan for AoU AMR
gw_data_aou_amr = get_gwas_nominal_hits("AoU_AMR",p_cutoff = 5e-8) # no hits nominal in NIAGADS
supp_fig1 = make_manhattan_easy_label(dataset = "AoU_AMR",NA,NA,fig_num="S1")

#######
# Supplemental Table 1: P values across cohorts for putatively causal variants.
supp_table_1 = get_p_put_causal_variants_across_studies() ; write_xlsx(supp_table_1,"Paper_Tables_Figures/supp_table_1.xlsx")

##########
# Figure 4a, i.e. Manhattan for AoU HISP
gw_data_aou_hisp = get_gwas_nominal_hits("AoU_HISP",p_cutoff = 1e-5) # no GW hits are nominal in NIA, so skip match
lead_variants_aou_hisp = get_lead_variants(gw_data_aou_hisp,get_gene_names=T)
figure4a = make_manhattan_easy_label(dataset = "AoU_HISP",gw_data_aou_hisp,
                                     lead_variants_aou_hisp,NA,
                                     fig_num="4A",show_gene_names=T,p_data_cut=1e-5,
                                     p_label_cut=1e-7)
loci = getNumberCommonRareLoci(aou_hisp_nominal_nia,maf_cut=0.01,return_lead=T)

#########
# Figure 4b, ie Manhattan for AoU AMR. All GW sig hits are not nominal in NIA
gw_data_aou_amr = get_gwas_nominal_hits("AoU_AMR",p_cutoff = 1e-5) 
lead_variants_aou_amr = get_lead_variants(gw_data_aou_amr,get_gene_names=T)
figure4b = make_manhattan_easy_label(dataset = "AoU_AMR",gw_data_aou_amr,
                                     lead_variants_aou_amr,NA,
                                     fig_num="4B",show_gene_names=T,p_data_cut = 1e-5,
                                     p_label_cut=1e-7)
loci = getNumberCommonRareLoci(aou_hisp_nominal_nia,maf_cut=0.01,return_lead=T)

##########
# Figure 4c, i.e. Manhattan for AoU NIA PROJ 3 SD. There is overlap with GW hits and NIAGADS nominal
# There are lots of loci here, so only show those loci that are nominal in NIAGADS
gw_data_aou_nia_proj = get_gwas_nominal_hits("AoU_NIA_PROJ_3SD",p_cutoff = 5e-8)
aou_nia_proj_nominal_nia = get_nia_nominal_hits(gw_data_aou_nia_proj)
lead_variants_aou_nia_proj = get_lead_variants(aou_nia_proj_nominal_nia,get_gene_names = F)
figure4c = make_manhattan_easy_label(dataset = "AoU_NIA_PROJ_3SD",gw_data_aou_nia_proj,
                                     lead_variants_aou_nia_proj,aou_nia_proj_nominal_nia,
                                     fig_num="4C",show_gene_names=T,p_data_cut = 5e-8,
                                     p_label_cut=5e-8)
loci = getNumberCommonRareLoci(aou_nia_proj_nominal_nia,maf_cut=0.01,return_lead=T)
