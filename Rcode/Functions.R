organize_gwas_data = function(pop) {
  if (pop %in% c("amr","eur","afr")) {
    files = list.files("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs",
                      pattern = glue("{pop}_multi_geno"),full.names = T)
    out_df = vroom(files[[1]],show_col_types = F) %>% filter(ID == "") # empty
    for (f in files) {
      if (!str_detect(f,"AD_any")) next
      print(glue("Reading {f}, file {which(files == f)} of {length(files)}.Focused on AD-by-proxy."))
      curr_f = vroom(f,show_col_types = F)
      out_df %<>% add_row(curr_f)
      gc()
    }
    vroom_write(out_df,glue("../Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_{pop}_summ_stats_AD_any.txt"))
  }else if (pop == "Proj_3SD") {
    files = list.files("../Data_Links/AoU_GWAS/AoU_NIAGADS_Proj_3SD",
                       pattern="regenie",full.names = T)
    out_df = vroom(files[[1]],show_col_types = F) %>% filter(ID == "") # empty
    for (f in files) {
      print(glue("Reading {f}, file {which(files == f)} of {length(files)}.Focused on AD-by-proxy."))
      curr_f = vroom(f,show_col_types = F) %>% mutate(Category = as.character(Category))
      out_df %<>% add_row(curr_f)
      gc()
    }
    vroom_write(out_df,glue("../Data_Links/AoU_GWAS/AoU_NIAGADS_Proj_3SD/aou_{pop}_summ_stats_AD_any.txt"))
  }else if (pop == 'MatchIt_Gen') {
    folder = "../Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_NoAgeSex_No_Anc_Covar"
    files = list.files(folder,
                       pattern="regenie",full.names = T)
    out_df = vroom(files[[1]],show_col_types = F) %>% filter(ID == "") # empty
    for (f in files) {
      print(glue("Reading {f}, file {which(files == f)} of {length(files)}.Focused on AD-by-proxy."))
      curr_f = vroom(f,show_col_types = F) %>% mutate(Category = as.character(Category))
      out_df %<>% add_row(curr_f)
      gc()
    }
    vroom_write(out_df %>% arrange(CHROM,GENPOS) %>% select(-Category),
                glue("{folder}/aou_{pop}_summ_stats_AD_any.txt"))
  }else if (pop == 'MatchIt_GenPhen') {
    folder = "../Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_AgeSex_No_Anc_Covar"
    files = list.files(folder,
                       pattern="regenie",full.names = T)
    out_df = vroom(files[[1]],show_col_types = F) %>% filter(ID == "") # empty
    for (f in files) {
      print(glue("Reading {f}, file {which(files == f)} of {length(files)}.Focused on AD-by-proxy."))
      curr_f = vroom(f,show_col_types = F) %>% mutate(Category = as.character(Category))
      out_df %<>% add_row(curr_f)
      gc()
    }
    vroom_write(out_df %>% arrange(CHROM,GENPOS) %>% select(-Category),
                glue("{folder}/aou_{pop}_summ_stats_AD_any.txt"))
  }
  print("Done aggregating statistics")
}
get_gwas_nominal_hits = function(dataset,p_cutoff=1e-5,window=500) {
  print("HWE data obtained similar to AoU UKB NIA paper, just using this study's IDs for P<1e-5")
  # FIRST, GET THE DATA FOR NOMINAL HITS
  if (dataset == "NIAGADS_HISP") {
    gw_df = vroom("/n/holystore01/LABS/tanzi_lab/Users/jwillett/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.p_5e-2.txt",
                  show_col_types = F) %>%
      rename(Pval = p,ID = marker,POS = pos,Freq1 = eaf,CHR=chr) %>% 
      mutate(CHR = str_replace(CHR,"chr",""),CHR=ifelse(CHR=="X","23",CHR),
             CHR = as.numeric(CHR)) %>% arrange(CHR,POS) %>%
      filter(Pval <= p_cutoff) %>%
      mutate(ID = glue("{CHR}-{POS}-{otherallele}-{effectallele}")) 
    hardy = vroom("working/anc_hwe_midp.txt") %>% 
      filter(MIDP_AFR > 1e-15,MIDP_AMR > 1e-15,MIDP_EAS > 1e-15, MIDP_EUR > 1e-15,
             MIDP_MID > 1e-15, MIDP_SAS > 1e-15) # remove AoU variants in HWE
    gw_df %<>% filter(ifelse(nchar(effectallele)==1 & nchar(otherallele)==1,
                             ID %in% hardy$ID | ID %in% hardy$IDrev,ID %in% hardy$ID))  %>%
      rename(Rsid=rsid,A1FREQ=Freq1,BETA=beta,SE=se)
  }else if (dataset == "AoU_NIA_Gen_GW") {
    gw_df = vroom("../Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_NoAgeSex_No_Anc_Covar/aou_MatchIt_Gen_summ_stats_AD_any_p_1e-2.txt",show_col_types = F) %>%
      mutate(Pval = 10^(-LOG10P)) %>% filter(Pval <= p_cutoff) %>% rename(CHR=CHROM,POS=GENPOS) %>%
      arrange(CHR,POS) %>% mutate(nearestgene = NA) %>%
      filter(A1FREQ*N >= 20,(1-A1FREQ)*N >= 20)
    hardy = vroom("working/anc_hwe_midp.txt") %>% 
      filter(MIDP_AFR > 1e-15,MIDP_AMR > 1e-15,MIDP_EAS > 1e-15, MIDP_EUR > 1e-15,
             MIDP_MID > 1e-15, MIDP_SAS > 1e-15) # remove AoU variants in HWE
    gw_df %<>% filter(ifelse(nchar(ALLELE0)==1 & nchar(ALLELE1)==1,
                             ID %in% hardy$ID | ID %in% hardy$IDrev,ID %in% hardy$ID))
  }else if (dataset == "AoU_NIA_GenPhen_GW") {
    gw_df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_AgeSex_No_Anc_Covar/aou_MatchIt_summ_stats_AD_any_p_1e-1.txt",show_col_types = F) %>%
      mutate(Pval = 10^(-LOG10P)) %>% filter(Pval <= p_cutoff) %>% rename(CHR=CHROM,POS=GENPOS) %>%
      arrange(CHR,POS) %>% mutate(nearestgene = NA) %>%
      filter(A1FREQ*N >= 20,(1-A1FREQ)*N >= 20)
    hardy = vroom("working/anc_hwe_midp.txt") %>% 
      filter(MIDP_AFR > 1e-15,MIDP_AMR > 1e-15,MIDP_EAS > 1e-15, MIDP_EUR > 1e-15,
             MIDP_MID > 1e-15, MIDP_SAS > 1e-15) # remove AoU variants in HWE
    gw_df %<>% filter(ifelse(nchar(ALLELE0)==1 & nchar(ALLELE1)==1,
                             ID %in% hardy$ID | ID %in% hardy$IDrev,ID %in% hardy$ID))
  }
  out_df = assignLociNumbers(gw_df %>% mutate(Locus=NA,NewOld = NA))
  
  # SECOND, CLASSIFY HITS AS NOVEL VS OLD (USING CODE FROM AOU VS UKB META ANALYSIS)
  known_loci = vroom("/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/known_loci_chrpos.txt")
  for (row in 1:nrow(out_df)) {
    curr_b_loc = known_loci %>% filter(CHR == out_df$CHR[[row]],
                                       POS - window*1000 <= out_df$POS[[row]],
                                       POS + window*1000 >= out_df$POS[[row]]) 
    if (nrow(curr_b_loc) > 0) out_df$NewOld[[row]] = "Old"
    else out_df$NewOld[[row]] = "New"
  }

  # FOURTH, add CHRPOS column to make intersecting datasets easier
  out_df %<>% mutate(CHRPOS = glue("{CHR}-{POS}"))
  
  return(out_df)
}
assignLociNumbers = function(df,cutoff=5e5) {
  # Number the loci
  curr_loc = 1
  out_df = df
  for (row in 1:nrow(df)) {
    r = out_df[row,]
    if (row == 1) { out_df$Locus[[row]] = curr_loc ; next }
    r_prev = out_df[(row-1),]
    if (is.na(r$Locus)) {
      if (r$CHR == r_prev$CHR & (r$POS-r_prev$POS <= cutoff)) { # joined with prior locus
        out_df$Locus[[row]] = r_prev$Locus
      }else{ # new locus
        curr_loc = curr_loc + 1
        out_df$Locus[[row]] = curr_loc
      }
    }
  }
  print(glue("Num total loci pre-QC: {length(unique(out_df$Locus))}"))
  return(out_df)
}
get_lead_variants = function(df,get_gene_names,meta = F) {
  out_df = df %>% filter(POS == -1) # so get empty data frame to add to 
  for (loc in unique(df$Locus)) {
    tmp = df %>% filter(Locus == loc)
    lead = tmp %>% filter(Pval == min(tmp$Pval))
    if (loc == 30) print(nrow(lead))
    out_df %<>% add_row(lead)
  }
  
  # When hits are not nominal in NIAGADS, still get gene names for plot
  if (get_gene_names) { # pull from FAVOR
    favor = vroom("../Data_Links/FAVOR/favor_hits_annot.txt",skip=1)
    names(favor) = names(vroom("../Data_Links/FAVOR/favor_hits_annot.txt",n_max = 0,delim='\t'))
    favor %<>% filter(variant_vcf %in% out_df$ID) %>% select(variant_vcf,genecode_comprehensive_info) %>%
      rename(ID=variant_vcf) %>% distinct(ID,.keep_all = T)
    
    # merge and deal with variants not in favor (gene pulled from VEP)
    merged = merge(out_df,favor,by='ID') %>% select(-nearestgene) %>%
      rename(nearestgene = genecode_comprehensive_info) %>%
      mutate(nearestgene = ifelse(ID == "17-75801070-ATTT-A",'UNK',nearestgene)) %>%
      mutate(nearestgene = ifelse(ID == "17-971686-GGTGA-G",'NXN',nearestgene)) %>%
      mutate(nearestgene = ifelse(ID == "19-22712186-CAT-C",'ZNF492-AS',nearestgene)) %>%
      mutate(nearestgene = ifelse(ID == "7-89881278-AT-A",'STEAP2-AS1',nearestgene)) %>%
      mutate(nearestgene = ifelse(ID == "12-121914347-C-CT",'PSMD9',nearestgene))
    
    # Print out any variants that did not make it
    print(glue("Num variants not in FAVOR: {nrow(merged %>% filter(ID %notin% out_df$ID))}"))
    out_df = merged
  }
  
  # Get rsids
  favor_rsids = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/FAVOR_Rsids.txt") %>%
    rename(ID = FAVOR_VCF)
  out_df_rsids = merge(out_df,favor_rsids,by='ID')
  print(glue('Nrow wo rsids: {nrow(out_df)}. Nrow w rsids: {nrow(out_df_rsids)}'))
  if (nrow(out_df) != nrow(out_df_rsids))
    print(glue("Missing rsids for: {(out_df %>% filter(ID %notin% out_df_rsids$ID))$ID}"))
  
  # flip major alleles
  for (row in 1:nrow(out_df_rsids)){
    curr_row = out_df_rsids[row,]
    if (!meta) {
      if (curr_row$A1FREQ >= 0.5) {
        out_df_rsids$BETA[[row]] = -1 * out_df_rsids$BETA[[row]]
        out_df_rsids$A1FREQ[[row]] = 1 - out_df_rsids$A1FREQ[[row]]
        a0 = out_df_rsids$ALLELE0[[row]] ; a1 = out_df_rsids$ALLELE0[[row]]
        out_df_rsids$ALLELE0[[row]] = a1
        out_df_rsids$ALLELE1[[row]] = a0
      }
    }else{
      if (curr_row$Freq1 >= 0.5) {
        out_df_rsids$Effect[[row]] = -1 * out_df_rsids$Effect[[row]]
        out_df_rsids$Freq1[[row]] = 1 - out_df_rsids$Freq1[[row]]
        a1 = out_df_rsids$Allele1[[row]] ; a2 = out_df_rsids$Allele2[[row]]
        out_df_rsids$Allele1[[row]] = a2
        out_df_rsids$Allele2[[row]] = a1
      }
    }
  }
  return(out_df_rsids %>% arrange(CHR,POS))
}
write_demographic_data_for_export = function() {
  demo_data_niagads = vroom("../Data_Links/NIAGADS_Personal/NIAGADS_demographic_data_all.txt") %>%
    filter(Ethnicity == "Hisp")
  out_df = demo_data_niagads %>% select(Age,Sex,superpop2)
  vroom_write(out_df,"figure_2_demographic_data.txt")
}
make_manhattan_easy_label = function(dataset,gw_qc,lead_var,fig_num,
                                     show_gene_names=T,only_new_old,
                                     p_data_cut,p_label_cut,show_highlights) {
  print("HWE QC already dealt with (incorporated per gw_qc filtering")
  if (dataset == "NIAGADS_HISP") {
    show_highlights = FALSE
    gw_df = vroom("/n/holystore01/LABS/tanzi_lab/Users/jwillett/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.p_5e-2.txt",
                  show_col_types = F) %>%
      rename(Pval = p,ID = marker,GENPOS = pos,Freq1 = eaf,CHR=chr) %>% 
      mutate(CHR = str_replace(CHR,"chr",""),CHR=ifelse(CHR=="X","23",CHR),
             CHR = as.numeric(CHR)) %>% arrange(CHR,GENPOS) %>%
      mutate(nearestgene = ifelse(ID %in% lead_var$ID,nearestgene,NA)) %>%
      mutate(Nominal = TRUE,NewOld = NA)
    for (var in 1:nrow(gw_qc)) {
      if (gw_qc$ID[[var]] %in% gw_df$ID)
        gw_df$ID[[which(gw_qc$ID[[var]] %in% gw_df$ID)]] = gw_qc$NewOld[[var]]
    }
  }else if (dataset == "AoU_NIA_Gen_GW") {
    gw_df = vroom("../Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_NoAgeSex_No_Anc_Covar/aou_MatchIt_Gen_summ_stats_AD_any_pval_p_1e-1.txt",show_col_types = F) %>%
      mutate(Pval = 10^(-LOG10P)) %>% mutate(nearestgene = NA,Nominal = FALSE,NewOld=NA) %>%
      rename(Freq1 = A1FREQ,CHR=CHROM) %>% arrange(CHR,GENPOS) %>%
      filter(ID %in% gw_qc$ID | Pval > p_label_cut)
    # add gene name and new old
    for (row in 1:nrow(lead_var)) {
      if (lead_var$Pval[[row]] > p_label_cut) next
      match = which(gw_df$ID == lead_var$ID[[row]])
      gw_df$nearestgene[[match]] = lead_var$nearestgene[[row]]
      gw_df$NewOld[[match]] = lead_var$NewOld[[row]]
    }
  }else if (dataset == "AoU_NIA_GenPhen_GW") {
    gw_df = vroom("../Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_AgeSex_No_Anc_Covar/aou_MatchIt_summ_stats_AD_any_p_1e-1.txt",show_col_types = F) %>%
      mutate(Pval = 10^(-LOG10P)) %>% mutate(nearestgene = NA,Nominal = FALSE,NewOld=NA) %>%
      rename(Freq1 = A1FREQ,CHR=CHROM) %>% arrange(CHR,GENPOS) %>%
      filter(ID %in% gw_qc$ID | Pval > p_label_cut)
    # add gene name and newold
    for (row in 1:nrow(lead_var)) {
      if (lead_var$Pval[[row]] > p_label_cut) next
      match = which(gw_df$ID == lead_var$ID[[row]])
      gw_df$nearestgene[[match]] = lead_var$nearestgene[[row]]
      gw_df$NewOld[[match]] = lead_var$NewOld[[row]]
    }
  }else if (dataset == "Meta_AoU_NIA_Gen_GW") {
    gw_df = vroom("meta/aou_nia_hisp_genmatch_metaanalysis_chrposrefalt_p_1e-1.TBL",show_col_types = F) %>%
      mutate(nearestgene = NA,Nominal = FALSE,CHR = as.numeric(CHR),NewOld=NA) %>%
      rename(GENPOS=POS,Pval=`P-value`,ID=MarkerName) %>% arrange(CHR,GENPOS) %>%
      filter(ID %in% gw_qc$ID | Pval > p_label_cut) 
    
    # add gene name and NewOld
    for (row in 1:nrow(lead_var)) {
      if (lead_var$Pval[[row]] > p_label_cut) next
      match = which(gw_df$ID == lead_var$ID[[row]])
      gw_df$nearestgene[[match]] = lead_var$nearestgene[[row]]
      gw_df$NewOld[[match]] = lead_var$NewOld[[row]]
    }
    
    # Focus on variants present in both studies
    print(glue("Sig variants wo HetDf filter: {nrow(gw_df %>% filter(Pval <= 5e-8))}"))
    gw_df %<>% filter(HetDf > 0)
    print(glue("Sig variants w HetDf filter: {nrow(gw_df %>% filter(Pval <= 5e-8))}"))
  }else if (dataset == "Meta_AoU_NIA_GenPhen_GW") {
    gw_df = vroom("meta/aou_nia_hisp_genphenmatch_metaanalysis_chrposrefalt_p_1e-1.TBL",show_col_types = F) %>%
      mutate(nearestgene = NA,Nominal = FALSE,CHR = as.numeric(CHR),NewOld=NA) %>%
      rename(GENPOS=POS,Pval=`P-value`,ID=MarkerName) %>% arrange(CHR,GENPOS) %>%
      filter(ID %in% gw_qc$ID | Pval > p_label_cut) 
    
    # add gene name
    for (row in 1:nrow(lead_var)) {
      if (lead_var$Pval[[row]] > p_label_cut) next
      match = which(gw_df$ID == lead_var$ID[[row]])
      gw_df$nearestgene[[match]] = lead_var$nearestgene[[row]]
      gw_df$NewOld[[match]] = lead_var$NewOld[[row]]
    }
    
    # Focus on variants present in both studies
    print(glue("Sig variants wo HetDf filter: {nrow(gw_df %>% filter(Pval <= 5e-8))}"))
    gw_df %<>% filter(HetDf > 0)
    print(glue("Sig variants w HetDf filter: {nrow(gw_df %>% filter(Pval <= 5e-8))}"))
  }
  
  # clean up gene name
  gw_df %<>% mutate(nearestgene = str_replace_all(nearestgene, "\\(dist.*?\\)", ""),
                    nearestgene = str_replace_all(nearestgene," ",""),
                    nearestgene = gsub("\\([^)]*\\)", "", nearestgene),
                    nearestgene = str_replace(nearestgene,'NONE,','')) 
  
  # show new or old
  if (only_new_old == "old") 
    gw_df %<>% mutate(nearestgene = ifelse(NewOld == "New","NL",nearestgene))
  else if (only_new_old == "new")
    gw_df %<>% mutate(nearestgene = ifelse(NewOld == "Old",NA,nearestgene))
  
  # all variants
  print('Making all variant plot')
  ggmanh::manhattan_plot(gw_df %>% filter(!Nominal) %>% add_row(gw_df %>% filter(Nominal)),
                         pval.colname = "Pval",chr.colname = "CHR",
                         pos.colname = "GENPOS",highlight.colname = "Nominal",
                         color.by.highlight=show_highlights,highlight.col = c("grey","red"),
                         label.colname = "nearestgene",rescale=T,
                         outfn = glue("Paper_Tables_Figures/Figure{fig_num}_all.png"),
                         max.overlaps = 20,label.font.size=2,signif = c(5e-08))
  
  # common variants
  print('Making common variant plots')
  ggmanh::manhattan_plot(gw_df %>% filter(Freq1 >= 0.01,Freq1 <= 0.99),
                         pval.colname = "Pval",chr.colname = "CHR",
                         pos.colname = "GENPOS",highlight.colname = "Nominal",
                         color.by.highlight = show_highlights,
                         label.colname = "nearestgene",rescale=T,
                         outfn = glue("Paper_Tables_Figures/Figure{fig_num}_common_0.01.png"),
                         max.overlaps = 1000,label.font.size=2,signif = c(5e-08))
  # 
  # # Rare variants
  print('Making rare variant plot')
  ggmanh::manhattan_plot(gw_df %>% filter(Freq1 < 0.01 | Freq1 > 0.99),
                         pval.colname = "Pval",chr.colname = "CHR",
                         pos.colname = "GENPOS",highlight.colname = "Nominal",
                         color.by.highlight = show_highlights,
                         label.colname = "nearestgene",rescale=T,
                         outfn = glue("Paper_Tables_Figures/Figure{fig_num}_rare_0.01.png"),
                         max.overlaps = 1000,label.font.size=2,signif = c(5e-08))
  
  print(glue("Made figure. Path: Figure{fig_num}_x.pdf"))
  return(gw_df)
}
intersect_for_p_values = function(lead_nia,lead_aou_genpcs,lead_aou_genpcs_phen) {
  chrpos = data.frame(CHRPOS = glue("{lead_nia$CHR}-{lead_nia$POS}")) %>%
    add_row(data.frame(CHRPOS = glue("{lead_aou_genpcs$CHR}-{lead_aou_genpcs$POS}"))) %>%
    add_row(data.frame(CHRPOS = glue("{lead_aou_genpcs_phen$CHR}-{lead_aou_genpcs_phen$POS}")))
  vroom_write(chrpos,"working/chrpos_for_gwas_intersect.txt")
  
  # then run the bash script for intersection
  system("bash intersect_gwas.sh")
}
make_gwas_lead_locus_table = function(lead_nia,lead_aou_genpcs,lead_aou_genpcs_phen) {
  out_df = data.frame(ID=NA,CHR=NA,Rsid=NA,Gene=NA,MAF=NA,EA=NA,OR=NA,CI=NA,NIA_P=NA,AoU_Gen_P=NA,
                      AoU_GenPhen_P=NA,NewOld=NA)
  for (lead in list(lead_nia,lead_aou_genpcs,lead_aou_genpcs_phen)) {
    out_df %<>% add_row(Rsid="")
    out_df %<>% add_row(ID=lead$ID,CHR=lead$CHR,Rsid=lead$Rsid,Gene=lead$nearestgene,
                        MAF=lead$A1FREQ,EA=lead$ALLELE1,OR=round(exp(lead$BETA),2),
                        CI=glue("({round(exp(lead$BETA-1.96*lead$SE),2)},{round(exp(lead$BETA+1.96*lead$SE),2)})"),
                        NIA_P=NA,AoU_Gen_P=NA,AoU_GenPhen_P=NA,NewOld=lead$NewOld)
  }
  # Next, get p-values for different studies.
  intersect_data = lapply(c("nia",'aou_genpconly','aou_genphen'),function(x) {
    vroom(glue("working/gwas_hits_{x}_intersect.txt"),show_col_types = F)
  })
  for (row in 1:nrow(out_df)) {
    if (is.na(out_df$ID[[row]])) next
    id = out_df$ID[[row]] ; id_split = str_split(id,"-")[[1]]
    idrev = glue("{id_split[[1]]}-{id_split[[2]]}-{id_split[[4]]}-{id_split[[3]]}")
    
    # NIAGADS
    nia_match = intersect_data[[1]] %>% mutate(marker = str_replace_all(marker,":","-")) %>%
      filter(marker == id | marker == idrev)
    if (nrow(nia_match)>0) out_df$NIA_P[[row]] = nia_match$p
    
    # AoU GEN PCs ONLY
    aou_genpcs = intersect_data[[2]] %>% filter(ID == id | ID == idrev)
    if (nrow(aou_genpcs)>0) out_df$AoU_Gen_P[[row]] = aou_genpcs$Pval
    
    # AoU GEN PCs PHEN
    aou_genpcs_phen = intersect_data[[3]] %>% filter(ID == id | ID == idrev)
    if (nrow(aou_genpcs_phen)>0) out_df$AoU_GenPhen_P[[row]] = aou_genpcs_phen$Pval
  }
  return(out_df)
}
make_meta_lead_locus_table = function(lead_genonly,lead_genphen) {
  out_df = data.frame(ID=NA,CHR=NA,Rsid=NA,Gene=NA,MAF=NA,EA=NA,OR=NA,CI=NA,MetaGenP=NA,MetaGenPhenP=NA,
                      NIA_P=NA,AoU_Gen_P=NA,
                      AoU_GenPhen_P=NA,NewOld=NA)
  for (lead in list(lead_genonly,lead_genphen)) {
    out_df %<>% add_row(Rsid="")
    out_df %<>% add_row(ID=lead$ID,CHR=as.numeric(lead$CHR),Rsid=lead$Rsid,Gene=lead$nearestgene,
                        MAF=lead$Freq1,EA=toupper(lead$Allele1),OR=round(exp(lead$Effect),2),
                        CI=glue("({round(exp(lead$Effect-1.96*lead$StdErr),2)},{round(exp(lead$Effect+1.96*lead$StdErr),2)})"),
                        MetaGenP=NA,MetaGenPhenP=NA,NIA_P=NA,AoU_Gen_P=NA,AoU_GenPhen_P=NA,
                        NewOld=lead$NewOld)
  }
  # Next, get p-values for different studies.
  intersect_data = lapply(c("meta_genonly","meta_genphen","gwas_nia","gwas_genonly",
                            "gwas_genphen"),function(x) {
    vroom(glue("working/meta_hits_{x}_intersect.txt"),show_col_types = F)
  })
  for (row in 1:nrow(out_df)) {
    if (is.na(out_df$ID[[row]])) next
    id = out_df$ID[[row]] ; id_split = str_split(id,"-")[[1]]
    idrev = glue("{id_split[[1]]}-{id_split[[2]]}-{id_split[[4]]}-{id_split[[3]]}")
    
    # GEN ONLY META CHECK
    meta_gen_match = intersect_data[[1]] %>% filter(MarkerName == id | MarkerName == idrev)
    if (nrow(meta_gen_match)>0) out_df$MetaGenP[[row]] = meta_gen_match$`P-value`
    
    # GEN PHEN META CHECK
    meta_genphen_match = intersect_data[[2]] %>% filter(MarkerName == id | MarkerName == idrev)
    if (nrow(meta_genphen_match)>0) out_df$MetaGenPhenP[[row]] = meta_genphen_match$`P-value`
    
    # NIAGADS
    nia_match = intersect_data[[3]] %>% mutate(marker = str_replace_all(marker,":","-")) %>%
      filter(marker == id | marker == idrev)
    if (nrow(nia_match)>0) out_df$NIA_P[[row]] = nia_match$p
    
    # AoU GEN PCs ONLY gwas
    aou_genpcs = intersect_data[[4]] %>% filter(ID == id | ID == idrev)
    if (nrow(aou_genpcs)>0) out_df$AoU_Gen_P[[row]] = aou_genpcs$Pval
    
    # AoU GEN PCs PHEN
    aou_genpcs_phen = intersect_data[[5]] %>% filter(ID == id | ID == idrev)
    if (nrow(aou_genpcs_phen)>0) out_df$AoU_GenPhen_P[[row]] = aou_genpcs_phen$Pval
  }
  return(out_df)
}


make_functional_table_niagads = function(df) {
  out_df = data.frame(ID=df$ID,CHR=df$CHR,rsid=df$rsid,Gene=NA,CADD_Phred=NA,
                      EpiActivMax=NA,EpiReprMax=NA,EpiTransMax=NA,Conserved=NA)
  
  # Get favor data
  favor_annot = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/favor_hits_annot.txt",skip=1,delim=",",col_names=F)
  names(favor_annot) = names(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/favor_hits_annot.txt",n_max=1,delim="\t"))
  
  # Match annot to out_df
  for (row in 1:nrow(out_df)) {
    if (row == 4) next
    favor_row = favor_annot %>% filter(variant_vcf == out_df$ID[[row]])
    if (nrow(favor_row)<1) next
    out_df$Gene[[row]] = favor_row$genecode_comprehensive_info
    out_df$CADD_Phred[[row]] = round(favor_row$cadd_phred,2)
    # out_df$EnhancerLinked[[row]] = ifelse(!is.na(favor_row$cage_enhancer) |  # none are enhancer linked
    #                                         !is.na(favor_row$genehancer) | 
    #                                         !is.na(favor_row$super_enhancer),
    #                                       "TRUE","FALSE")
    out_df$EpiActivMax[[row]] = round(max(c(favor_row$apc_epigenetics_active,favor_row$encodeh3k27ac_sum,
                                      favor_row$encodeh3k4me1_sum,favor_row$encodeh3k4me2_sum,
                                      favor_row$encodeh3k4me3_sum,favor_row$encodeh3k9ac_sum,
                                      favor_row$encodeh4k20me1_sum,favor_row$encodeh2afz_sum),na.rm=T),2)
    out_df$EpiReprMax[[row]] = round(max(c(favor_row$apc_epigenetics_repressed,favor_row$encodeh3k9me3_sum,
                                     favor_row$encodeh3k27me3_sum),na.rm=T),2)
    out_df$EpiTransMax[[row]] = round(max(c(favor_row$apc_epigenetics_transcription,
                                      favor_row$encodeh3k36me3_sum,
                                      favor_row$encodeh3k79me2_sum),na.rm=T),2)
    out_df$Conserved[[row]] = ifelse((!is.na(favor_row$priphylop) & favor_row$priphylop >= 0.3) |
                                       (!is.na(favor_row$priphcons) & favor_row$priphcons >= 0.3) |
                                       (!is.na(favor_row$mamphcons) & favor_row$mamphcons >= 0.3) |
                                       (!is.na(favor_row$verphcons) & favor_row$verphcons >= 0.3) |
                                       (!is.na(favor_row$gerp_n) & favor_row$gerp_n >= 10) |
                                       (!is.na(favor_row$gerp_s) & favor_row$gerp_s >= 10),
                                     "TRUE","FALSE")
  }
  return(out_df)
}
getSingleCellData = function(p_cutoff) {
  # results wrt comparing larger groups, ie stratified by pathological evidence vs cognitive impairment
  big_picture_results = vroom("../00_AoU/scrna_logfold_stats_by_cell_group_by_pathology_and_symptoms_1234.txt") %>%
    filter(PValAdj <= p_cutoff) %>% arrange(Gene,desc(AvgLog2FC)) %>% group_by(Gene) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"second_vs_first"),"2v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"third_vs_first"),"3v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_first"),"4v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"third_vs_second"),"3v2",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_second"),"4v2",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_third"),"4v3",Comparison)) 
  
  # results wrt cognitive impairment in the setting of different pathological features
  small_picture_results = vroom("../00_AoU/scrna_logfold_stats_by_cell.txt") %>%
    filter(PValAdj <= p_cutoff) %>% arrange(Gene,desc(AvgLog2FC)) %>% 
    mutate(Comparison = str_replace(Comparison,".csv","")) %>% group_by(Gene) %>%
    mutate(CellPop = ifelse(str_detect(Comparison,"second_vs_first"),glue("{CellPop} 2v1"),CellPop)) %>%
    mutate(CellPop = ifelse(str_detect(Comparison,"third_vs_second"),glue("{CellPop} 3v2"),CellPop))
  
  out_list = list(big_picture_results,small_picture_results)
  names(out_list) = c("Gross_Phenotypes","Micro_Phenotypes")
  return(out_list)
}
make_general_expr_table = function(df_nia,df_aou,sc_data) {
  # Only care about the genes, but need variants to keep track of epimax
  df_nia_processed = df_nia %>% filter(EpiActivMax >= 20 | EpiReprMax >= 20 | EpiTransMax >= 20)
  df_aou_processed = df_aou %>% filter(EpiActivMax >= 20 | EpiReprMax >= 20 | EpiTransMax >= 20 | !is.na(EnhancerGene)) %>%
    arrange(CHR)
  out_df = data.frame(ID = df_nia_processed$ID,CHR=df_nia_processed$CHR,rsid=df_nia_processed$rsid,
                      Gene=df_nia_processed$Gene,EnhancerGene=NA,MaxLocusEpiPhred=NA,CellwComp=NA) %>% 
    add_row() %>%
    add_row(data.frame(ID = df_aou_processed$ID,CHR=df_aou_processed$CHR,rsid=df_aou_processed$rsid,
                       Gene=df_aou_processed$Gene,EnhancerGene=df_aou_processed$EnhancerGene,
                       MaxLocusEpiPhred=NA,CellwComp=NA))

  for (row in 1:nrow(out_df)) {
    if (is.na(out_df$ID[[row]])) next
    # First get the genes on the row
    row_genes = out_df[row,] %>% select(Gene,EnhancerGene) %>% 
      mutate(Gene = glue("{Gene},{EnhancerGene}")) %>%
      separate_rows(Gene, sep = ",") %>% filter(!is.na(Gene),Gene != 'NA')
    
    # Then clean up gene names to enable matching
    row_genes %<>% mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
                       Gene = str_replace_all(Gene," ",""),
                       Gene = gsub("\\([^)]*\\)", "", Gene),
                       Gene = str_replace(Gene,'NONE,',''))

    # then get the expr change
    genes_with_expr = character()
    for (g in unique(row_genes$Gene)) { # cross reference gene, unique to avoid duplicates
      sc_res_match = sc_data %>% filter(Gene == g) %>% mutate(CellPopComp = paste(CellPop,Comparison))
      if (nrow(sc_res_match)>0) {
        genes_with_expr %<>% append(g)
        out_df$CellwComp[[row]] = 
          paste0(out_df$CellwComp[[row]],paste(glue("({g})"),
                                               glue("({sc_res_match$CellPopComp};{sc_res_match$PValAdj};{sc_res_match$AvgLog2FC})"),
                                                            collapse=","),collapse=",")
        if (length(unique(row_genes$Gene)) > 1 & 
            g != unique(row_genes$Gene)[[length(unique(row_genes$Gene))]])
          out_df$CellwComp[[row]] = paste0(out_df$CellwComp[[row]],",")
      }
    }
    # clean up gene expr change entry
    out_df$CellwComp[[row]] = clean_cell_expr_column(out_df[[row,"CellwComp"]],row_genes$Gene)
    if (length(genes_with_expr) == 1) {
      out_df$CellwComp[[row]] = str_replace_all(out_df$CellwComp[[row]],glue("\\({genes_with_expr[[1]]}\\) "),"")
    }
    
    # Then add the epi phred max value for further information on links
    if (out_df$ID[[row]] %in% df_nia_processed$ID) {
      match = df_nia_processed %>% filter(ID == out_df$ID[[row]])
      out_df$MaxLocusEpiPhred[[row]] = max(c(match$EpiActivMax,match$EpiReprMax,
                                             match$EpiTransMax))
    }else{
      match = df_aou_processed %>% filter(ID == out_df$ID[[row]])
      out_df$MaxLocusEpiPhred[[row]] = max(c(match$EpiActivMax,match$EpiReprMax,
                                             match$EpiTransMax))
    }
  }
  return(out_df)
}
clean_cell_expr_column = function(cell_entry,genes) { # genes are genes in the Gene column
  cell_entry = str_replace(cell_entry,"NA\\(","\\(")
  clean_str = ""

  spl = str_split(cell_entry,",")[[1]] # get each separate entry
  for (gene in genes) { # to deal with more than one gene
    gene_edit = str_replace_all(gene,"\\(dist.*?\\)", "")
    gene_edit = str_replace_all(gene_edit," ","")
    gene_edit = gsub("\\([^)]*\\)", "", gene_edit)
    
    spl_g = grep(glue("\\({gene_edit}\\)"),spl,value=T) # focus on entries for given gene
    
    cells = c("Exc ","Inh ","Oli ","OPC ","Ast ","Mic ")
    comps = c(" 2v1"," 3v1"," 4v1"," 3v2"," 4v2"," 4v3")
    
    for (cell in cells) {
      for (comp in comps) {
        count = length(which(str_detect(spl_g,cell) & str_detect(spl_g,comp)))
        # Also want to pull by p-value and the log2fc (2nd and 3rd entry delim by semicolons)
        if (count > 0) {
          # Sign sometimes vary by population, so less informative, thus include only when selected
          if (TRUE) {
            sign = paste0(get_logfc_sign(spl_g,cell,comp))
          }else sign = ""
          clean_str = glue("{clean_str} \\({gene_edit}\\) {cell} {comp} (x{count}) {sign}")
        }
      }
    }
    
    out_cell = str_squish(clean_str)
    out_cell = gsub("\\\\","",out_cell)
    
    if (length(genes) == 1) {
      out_cell = str_replace_all(out_cell,glue("\\({genes[[1]]}\\) "),"")
    }
  }
  return(out_cell)
}
get_logfc_sign = function(v,ce,co) { # take in vector of strings, cell, and comparison
  lines = v[which(str_detect(v,ce) & str_detect(v,co))]
  if (length(lines) < 1) return("")
  out_str = ""
  for (line in lines) {
    spl_l = str_split(line,";")[[1]][[3]]
    if (str_detect(spl_l,"-")) out_str = paste0(out_str,"-")
    else out_str = paste0(out_str,"+")
  }
  return(out_str)
}
make_path_specific_expr_table = function(df,sc_data) {
  out_df = data.frame(ID=NA,Path=NA,CHR=NA,rsid=NA,ProximalGene=NA,
                      MaxLocusEpiPhred=NA,CellwComp=NA)
  
  for (row in 1:nrow(df)) {
    if (row > 3) next
    curr_row = df[row,]
    
    # First get the genes on the row
    row_genes = curr_row %>% select(Gene) %>%
      mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
             Gene = str_replace_all(Gene," ",""),
             Gene = gsub("\\([^)]*\\)", "", Gene)) %>%
      separate_rows(Gene, sep = ",")
    
    # then get the expr change
    for (g in unique(row_genes$Gene)) { # cross reference gene, unique to avoid duplicates
      amyloid = gpath = nft = plaqd = plaqn = tangles = ""
      sc_res_match = sc_data %>% filter(Gene == g) 
      amyloid_match = sc_res_match %>% filter(Outcome == "amyloid") 
      gpath_match = sc_res_match %>% filter(Outcome == "gpath") 
      nft_match = sc_res_match %>% filter(Outcome == "nft") 
      plaqd_match = sc_res_match %>% filter(Outcome == "plaq_d")
      plaqn_match = sc_res_match %>% filter(Outcome == "plaq_n")
      tangles_match = sc_res_match %>% filter(Outcome == "tangles")
      
      if (nrow(amyloid_match)>0) amyloid = paste(amyloid,paste(glue("({g})"),
        glue("({amyloid_match$CellPop};{amyloid_match$PValAdj};{amyloid_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(gpath_match)>0) gpath = paste(gpath,paste(glue("({g})"),
        glue("({gpath_match$CellPop};{gpath_match$PValAdj};{gpath_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(nft_match)>0) nft = paste(nft,paste(glue("({g})"),
        glue("({nft_match$CellPop};{nft_match$PValAdj};{nft_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(plaqd_match)>0) plaqd = paste(plaqd,paste(glue("({g})"),
        glue("({plaqd_match$CellPop};{plaqd_match$PValAdj};{plaqd_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(plaqn_match)>0) plaqn = paste(plaqn,paste(glue("({g})"),
        glue("({plaqn_match$CellPop};{plaqn_match$PValAdj};{plaqn_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(tangles_match)>0) tangles = paste(tangles,paste(glue("({g})"),
        glue("({tangles_match$CellPop};{tangles_match$PValAdj};{tangles_match$AvgLog2FC})"),collapse=","),",")
      
      if (amyloid != "") {
        amyloid_out = clean_cell_expr_column(amyloid,g)
        out_df %<>% add_row(data.frame(ID=curr_row$ID,Path="Amyloid",CHR=curr_row$CHR,
                                       rsid=curr_row$rsid,ProximalGene=g,
                                       MaxLocusEpiPhred=NA,CellwComp=amyloid_out))
      }
      if (gpath != "") {
        gpath_out = clean_cell_expr_column(gpath,g)
        out_df %<>% add_row(data.frame(ID=curr_row$ID,Path="Gpath",CHR=curr_row$CHR,
                                       rsid=curr_row$rsid,ProximalGene=g,
                                       MaxLocusEpiPhred=NA,CellwComp=gpath_out))
      }
      if (nft != "") {
        nft_out = clean_cell_expr_column(nft,g)
        out_df %<>% add_row(data.frame(ID=curr_row$ID,Path="NFT",CHR=curr_row$CHR,
                                       rsid=curr_row$rsid,ProximalGene=g,
                                       MaxLocusEpiPhred=NA,CellwComp=nft_out))
      }
      if (plaqd != "") {
        plaqd_out = clean_cell_expr_column(plaqd,g)
        out_df %<>% add_row(data.frame(ID=curr_row$ID,Path="PlaqD",CHR=curr_row$CHR,
                                       rsid=curr_row$rsid,ProximalGene=g,
                                       MaxLocusEpiPhred=NA,CellwComp=plaqd_out))
      }
      if (plaqn != "") {
        plaqn_out = clean_cell_expr_column(plaqn,g)
        out_df %<>% add_row(data.frame(ID=curr_row$ID,Path="PlaqN",CHR=curr_row$CHR,
                                       rsid=curr_row$rsid,ProximalGene=g,
                                       MaxLocusEpiPhred=NA,CellwComp=plaqn_out))
      }
      if (tangles != "") {
        tangles_out = clean_cell_expr_column(tangles,g)
        out_df %<>% add_row(data.frame(ID=curr_row$ID,Path="Tangles",CHR=curr_row$CHR,
                                       rsid=curr_row$rsid,ProximalGene=g,
                                       MaxLocusEpiPhred=NA,CellwComp=tangles_out))
      }
      out_df$MaxLocusEpiPhred[which(out_df$ProximalGene == g)] = 
        max(c(curr_row$EpiActivMax,curr_row$EpiReprMax,
              curr_row$EpiTransMax))
    }
  }
  
  
  return(out_df %>% drop_na(ID) %>% arrange(CHR) %>% select(-CHR,-ID))
}
get_variants_in_ld = function(dataset,r2_cutoff = 0.1) {
  # first get all the IDs
  put_causal_variants_nia = c("18:11144178:C:T","18:11144561:A:T",
                          "18:11144647:G:A","18:11145132:A:T")
  put_causal_variants_aou = str_replace_all(put_causal_variants_nia,":","-")
  if (dataset == "NIAGADS") { # using NIAGADS (smaller cohort) for reference LD
    ld_data = vroom("piezo2_r2_phased_matrix.vcor") %>%
      filter(ID_A %in% put_causal_variants_nia | ID_B %in% put_causal_variants_nia,PHASED_R2 >= r2_cutoff)
    variants_to_evaluate = unique(put_causal_variants_nia %>% append(ld_data$ID_A) %>% append(ld_data$ID_B))
    summ_stats_piezo2 = vroom("../Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.piezo2.txt",show_col_types = F)
    ordered_stats = summ_stats_piezo2 %>% filter(marker %in% put_causal_variants_nia) %>%
      add_row() %>% # empty line
      add_row(summ_stats_piezo2 %>% filter(marker %notin% put_causal_variants_nia,marker %in% variants_to_evaluate)) %>%
      rename(ID = marker,Pval = p)
  }else if (str_detect(dataset,"AoU")) {
    if (dataset == "AoU_HISP") {
      ld_data = vroom("aou_ld_data_piezo2/aou_hisp_piezo2_r2_1e-1.vcor",show_col_types = F)
      summ_stats_piezo2 = vroom("../Data_Links/AoU_GWAS/HISP_GWAS_Hisp_PCs/aou_hisp_hisp_pcs_geno_1e-1_mac_20_piezo2.txt",show_col_types = F)
    }else if (dataset == "AoU_AMR") {
      ld_data = vroom("aou_ld_data_piezo2/aou_amr_piezo2_r2_1e-1.vcor",show_col_types = F)
      summ_stats_piezo2 = vroom("../Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_amr_summ_stats_AD_any_piezo2.txt",show_col_types = F)
    }else if (dataset == "AoU_NIA_PROJ") {
      ld_data = vroom("aou_ld_data_piezo2/aou_nia_proj_3sd_piezo2_r2_1e-1.vcor",show_col_types = F)
      summ_stats_piezo2 = vroom("../Data_Links/AoU_GWAS/AoU_NIAGADS_Proj_3SD/aou_Proj_3SD_summ_stats_AD_any_piezo2.txt",show_col_types = F)
    }else if (dataset == "AoU_All") {
      ld_data = vroom("aou_ld_data_piezo2/aou_all_piezo2_r2_1e-1.vcor",show_col_types = F)
      summ_stats_piezo2 = vroom("../Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_piezo2locus.txt",show_col_types = F)
    }
    ld_data %<>% filter(ID_A %in% put_causal_variants_aou | ID_B %in% put_causal_variants_aou,PHASED_R2 >= r2_cutoff)
    variants_to_evaluate = unique(put_causal_variants_aou %>% append(ld_data$ID_A) %>% append(ld_data$ID_B))
    ordered_stats = summ_stats_piezo2 %>% filter(ID %in% put_causal_variants_aou) %>%
      add_row() %>% # empty line
      add_row(summ_stats_piezo2 %>% filter(ID %notin% put_causal_variants_aou,ID %in% variants_to_evaluate)) %>%
      mutate(Pval = 10^(-LOG10P)) 
  }else if (dataset == 'UKB') {
    ld_data = vroom("aou_ld_data_piezo2/ukb_all_piezo2_geno_1e-1_r2_1e-1.vcor",show_col_types = F)
    summ_stats_piezo2 = vroom("../Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_piezo2.regenie",show_col_types = F)
    ld_data %<>% filter(ID_A %in% put_causal_variants_aou | ID_B %in% put_causal_variants_aou,PHASED_R2 >= r2_cutoff)
    variants_to_evaluate = unique(put_causal_variants_aou %>% append(ld_data$ID_A) %>% append(ld_data$ID_B))
    ordered_stats = summ_stats_piezo2 %>% filter(ID %in% put_causal_variants_aou) %>%
      add_row() %>% # empty line
      add_row(summ_stats_piezo2 %>% filter(ID %notin% put_causal_variants_aou,ID %in% variants_to_evaluate)) %>%
      mutate(Pval = 10^(-LOG10P)) %>% rename(CHROM = `#CHROM`)
  }

  # return a data frame with the necessary information
  out_df = data.frame(ID=ordered_stats$ID,Rsid=NA,Pval = ordered_stats$Pval)

  return(out_df %>% mutate(ID = str_replace_all(out_df$ID,":","-")))
}
check_hits_ld_hits_all_datasets = function() {
  cohorts = c("NIAGADS","AoU_HISP","AoU_AMR","AoU_NIA_PROJ","AoU_All","UKB")
  ld_data = list() ; all_ids = character()
  for (cohort in cohorts) {
    print(glue("Getting LD data for {cohort}"))
    ld_data[[cohort]] = get_variants_in_ld(cohort)
    all_ids %<>% append(ld_data[[cohort]]$ID)
  }
  
  # now organize a common data frame
  out_df = data.frame(ID = unique(all_ids),NIAGADS_P=NA,AoU_HISP_P=NA,AoU_AMR_P=NA,
                      AoU_NIA_PROJ_P=NA,AoU_All_P=NA,UKB_P=NA)
  for (cohort in cohorts) {
    print(glue("On cohort: {cohort}"))
    if (cohort=="NIAGADS") summstats = vroom("../Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.piezo2.txt",show_col_types = F) %>%
        rename(CHROM=chr,GENPOS=pos,ALLELE0=otherallele,ALLELE1=effectallele) %>%
        mutate(LOG10P = -log10(p))
    else if (cohort == "AoU_HISP") summstats = vroom("../Data_Links/AoU_GWAS/HISP_GWAS_Hisp_PCs/aou_hisp_hisp_pcs_geno_1e-1_mac_20_piezo2.txt",show_col_types = F)
    else if (cohort == "AoU_AMR") summstats = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_amr_summ_stats_AD_any_piezo2.txt",show_col_types = F)
    else if (cohort == "AoU_NIA_PROJ") summstats = vroom("../Data_Links/AoU_GWAS/AoU_NIAGADS_Proj_3SD/aou_Proj_3SD_summ_stats_AD_any_piezo2.txt",show_col_types = F)
    else if (cohort == "AoU_All") summstats = vroom("../Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_piezo2locus.txt",show_col_types = F)
    else if (cohort == "UKB") summstats = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_piezo2.regenie",show_col_types = F) %>%
        rename(CHROM = `#CHROM`)
    
    summstats %<>% mutate(ID = glue("{CHROM}-{GENPOS}-{ALLELE0}-{ALLELE1}"), IDrev = glue("{CHROM}-{GENPOS}-{ALLELE1}-{ALLELE0}"))
    
    for (row in 1:nrow(out_df)) {
      curr_row = out_df[row,] ; curr_col = which(str_detect(names(out_df),cohort))
      if (curr_row$ID %in% summstats$ID | curr_row$ID %in% summstats$IDrev)
        rows = which(summstats$ID %in% curr_row$ID | summstats$IDrev %in% curr_row$ID)
        if (length(rows)>1) rows = which(summstats$ID %in% curr_row$ID)
        out_df[[curr_col]][[row]] = 10^(-1*summstats$LOG10P[rows])
    }
  }
  
  # isolate output with nominal significance in one of the other datasets
  return_df = out_df %>% filter(AoU_HISP_P <= 0.05 & (AoU_AMR_P <= 0.05 | 
                                AoU_NIA_PROJ_P <= 0.05 | AoU_All_P <= 0.05 | UKB_P <= 0.05))
  return(return_df)
}
get_nia_nominal_hits = function(df) {
  print('Reading in NIAGADS summary stats')
  nia_hisp_nominal_gwas = vroom("../Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.p_5e-2.txt",show_col_types = F) %>%
    mutate(CHRPOS=glue("{chr}-{pos}"),ID=glue("{chr}-{pos}-{otherallele}-{effectallele}"),
           IDrev=glue("{chr}-{pos}-{effectallele}-{otherallele}"))
  
  # get the nominal hits
  print('Filtering AoU hits by nominal')
  out_df = df %>% mutate(CHRPOS=glue("{CHR}-{POS}")) %>% filter(CHRPOS %in% nia_hisp_nominal_gwas$CHRPOS) %>%
    filter(ID %in% nia_hisp_nominal_gwas$ID | ID %in% nia_hisp_nominal_gwas$IDrev) %>%
    mutate(Rsid=NA,nearestgene=NA)
  if (nrow(out_df) == 0) stop("No variants in this dataset are nominal in NIA HISP")
  
  # add the rsids and gene names, where available
  for (row in 1:nrow(out_df)) {
    nia_match = nia_hisp_nominal_gwas %>% filter(CHRPOS == out_df$CHRPOS[[row]])
    out_df$nearestgene[[row]] = nia_match$nearestgene
    if (nia_match$ID == out_df$ID[[row]])
      out_df$Rsid[[row]] = nia_match$rsid
  }
  return(out_df)
}
getNumberCommonRareLoci = function(df,maf_cut = 0.005,return_lead) {
  num_common = num_rare = 0
  num_new_common = num_new_rare = 0
  
  out_df = df %>% filter(Locus == -1) # so get empty data frame to add to 
  for (loc in unique(df$Locus)) {
    tmp = df %>% filter(Locus == loc) # all are GW significant
    
    if (return_lead) {
      lead = tmp %>% filter(Pval == min(tmp$Pval))
      out_df %<>% add_row(lead)
      if (max(lead$Freq1) >= maf_cut) { # testing lead variant
        num_common = num_common + 1
        if (str_detect(lead$NewOld,"New")) num_new_common = num_new_common + 1
      }else {
        num_rare = num_rare + 1
        if (str_detect(lead$NewOld,"New")) num_new_rare = num_new_rare + 1
      }
    }else{
      if (nrow(tmp %>% filter(Freq1 >= maf_cut)) > 0) { # testing all variants in locus
        num_common = num_common + 1
        if (TRUE %in% str_detect(tmp$NewOld,"New")) num_new_common = num_new_common + 1
      }else {
        num_rare = num_rare + 1
        if (TRUE %in% str_detect(tmp$NewOld,"New")) num_new_rare = num_new_rare + 1
      }
    }
  }
  print(glue("Total num variants: {nrow(df)}"))
  print(glue("Num common loci: {num_common}. Num rare loci: {num_rare}"))
  print(glue("Num new common loci: {num_new_common}. Num new rare loci: {num_new_rare}"))
  print("Common locus: has at least one common variant in the locus, which could not be the lead")
  print("Rare locus: has ALL rare variants in the locus, including the lead")
  return(out_df)
}
get_aou_cohort_intersects = function(chrpos) {
  intersect_list = list()
  cohorts = c("AoU_HISP","AoU_AMR","AoU_NIA_PROJ_3SD",'AoU_NIA_MatchIt',"NIAGADS")
  intersect_list = lapply(cohorts,function(c) {
    print(glue("Curr cohort: {c}"))
    gc()
    if (c == "AoU_HISP")
      vroom("../Data_Links/AoU_GWAS/HISP_GWAS_Hisp_Pcs_Anc_Covar/aou_hisp_anc_covar.txt",show_col_types = F) %>%
      mutate(CHRPOS = glue("{CHROM}-{GENPOS}")) %>% filter(CHRPOS %in% chrpos)
    else if (c == "AoU_AMR")
      vroom("../Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_amr_summ_stats_AD_any.txt",show_col_types = F) %>%
      mutate(CHRPOS = glue("{CHROM}-{GENPOS}")) %>% filter(CHRPOS %in% chrpos)
    else if (c == "AoU_NIA_PROJ_3SD")
      vroom("../Data_Links/AoU_GWAS/AoU_NIAGADS_Proj_3SD/aou_Proj_3SD_summ_stats_AD_any.txt",show_col_types = F) %>%
      mutate(CHRPOS = glue("{CHROM}-{GENPOS}")) %>% filter(CHRPOS %in% chrpos)
    else if (c == "AoU_NIA_MatchIt")
      vroom("../Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_No_Anc_Covar/aou_MatchIt_summ_stats_AD_any.txt",show_col_types = F) %>%
      mutate(CHRPOS = glue("{CHROM}-{GENPOS}")) %>% filter(CHRPOS %in% chrpos)
    else if (c == "NIAGADS")
      vroom("../Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.dn8.gz",show_col_types = F) %>%
      mutate(CHRPOS = glue("{chr}-{pos}")) %>% filter(CHRPOS %in% chrpos)
  })
  names(intersect_list) = cohorts
  return(intersect_list)
}
make_aou_lead_variant_table = function(df_hisp,df_amr,df_proj,df_matchit,intersects) {
  # First organize input data
  out_df = data.frame(ID=df_hisp$ID,BETA=df_hisp$BETA,SE=df_hisp$SE,CHR=df_hisp$CHR,
                      Rsid=df_hisp$Rsid,ProximalGene=df_hisp$nearestgene,
                      MAF=round(df_hisp$Freq1,3),OR=round(exp(df_hisp$BETA),2),
                      CI=glue("({round(exp(df_hisp$BETA-1.96*df_hisp$SE),2)},{round(exp(df_hisp$BETA+1.96*df_hisp$SE),2)})"),
                      NewOld=df_hisp$NewOld,AoU_HISP_P=df_hisp$Pval,AoU_AMR_P=NA,
                      AoU_NIA_Similar_P = NA,AoU_NIA_MatchIt_P=NA,NIA_P=NA) %>% add_row() %>%
    add_row(data.frame(ID=df_amr$ID,BETA=df_amr$BETA,SE=df_amr$SE,CHR=df_amr$CHR,
                       Rsid=df_amr$Rsid,ProximalGene=df_amr$nearestgene,
                       MAF=round(df_amr$Freq1,3),OR=round(exp(df_amr$BETA),2),
                       CI=glue("({round(exp(df_amr$BETA-1.96*df_amr$SE),2)},{round(exp(df_amr$BETA+1.96*df_amr$SE),2)})"),
                       NewOld=df_amr$NewOld,AoU_HISP_P=NA,AoU_AMR_P=df_amr$Pval,
                       AoU_NIA_Similar_P = NA,AoU_NIA_MatchIt_P=NA,NIA_P=NA)) %>% add_row() %>%
    add_row(data.frame(ID=df_proj$ID,BETA=df_proj$BETA,SE=df_proj$SE,CHR=df_proj$CHR,
                       Rsid=df_proj$Rsid,ProximalGene=df_proj$nearestgene,
                       MAF=round(df_proj$Freq1,3),OR=round(exp(df_proj$BETA),2),
                       CI=glue("({round(exp(df_proj$BETA-1.96*df_proj$SE),2)},{round(exp(df_proj$BETA+1.96*df_proj$SE),2)})"),
                       NewOld=df_proj$NewOld,AoU_HISP_P=NA,AoU_AMR_P=NA,
                       AoU_NIA_Similar_P = df_proj$Pval,AoU_NIA_MatchIt_P = NA,NIA_P=NA)) %>% add_row() %>%
     add_row(data.frame(ID=df_matchit$ID,BETA=df_matchit$BETA,SE=df_matchit$SE,CHR=df_matchit$CHR,
                        Rsid=df_matchit$Rsid,ProximalGene=df_matchit$nearestgene,
                       MAF=round(df_matchit$Freq1,3),OR=round(exp(df_matchit$BETA),2),
                       CI=glue("({round(exp(df_matchit$BETA-1.96*df_matchit$SE),2)},{round(exp(df_matchit$BETA+1.96*df_matchit$SE),2)})"),
                       NewOld=df_matchit$NewOld,AoU_HISP_P=NA,AoU_AMR_P=NA,
                       AoU_NIA_Similar_P = NA,AoU_NIA_MatchIt_P=df_matchit$Pval,NIA_P=NA))
  
  # then fill in missing p values
  # PUT IN INTERSECT CODE
  for (row in 1:nrow(out_df)) {
    curr_row = out_df[row,] ; spl = str_split(curr_row$ID,"-")[[1]]
    if (is.na(curr_row$ID)) next
    id_rev = glue("{spl[[1]]}-{spl[[2]]}-{spl[[4]]}-{spl[[3]]}")
    
    if (is.na(curr_row$AoU_HISP_P)) out_df$AoU_HISP_P[[row]] = 10^(-1*(intersects$AoU_HISP %>% filter(ID == curr_row$ID))$LOG10P)
    if (is.na(curr_row$AoU_AMR_P)) {
      match = (intersects$AoU_AMR %>% filter(ID == curr_row$ID))
      if (nrow(match)>0) out_df$AoU_AMR_P[[row]] = 10^(-1*match$LOG10P)
    }
    if (is.na(curr_row$AoU_NIA_Similar_P)) out_df$AoU_NIA_Similar_P[[row]] = 10^(-1*(intersects$AoU_NIA_PROJ_3SD %>% filter(ID == curr_row$ID))$LOG10P)
    if (is.na(curr_row$AoU_NIA_MatchIt_P)) out_df$AoU_NIA_MatchIt_P[[row]] = 10^(-1*(intersects$AoU_NIA_MatchIt %>% filter(ID == curr_row$ID))$LOG10P)
    if (is.na(curr_row$NIA_P)) {
      nia_match = intersects$NIAGADS %>% mutate(marker = str_replace_all(marker,":","-")) %>% 
        filter(marker == curr_row$ID | marker == id_rev)
      if (nrow(nia_match)>0) out_df$NIA_P[[row]] = nia_match$p
    }
  }
  # Clean up gene names
  out_df %<>% mutate(ProximalGene = str_replace_all(ProximalGene, "\\(dist.*?\\)", ""),
                     ProximalGene = str_replace_all(ProximalGene," ",""),
                     ProximalGene = gsub("\\([^)]*\\)", "", ProximalGene),
                     ProximalGene = str_replace(ProximalGene,'NONE,','')) 
  
  # Flip major to minor if present
  # Flip alleles if not minor
  for (row in 1:nrow(out_df)) {
    curr_row = out_df[row,]
    if (is.na(curr_row$MAF)) next
    if (curr_row$MAF > 0.5) {
      out_df$MAF[[row]] = 1 - curr_row$MAF
      curr_row$BETA = curr_row$BETA * -1
      out_df$OR[[row]] = round(exp(curr_row$BETA),2)
      out_df$CI[[row]] = glue("({round(exp(curr_row$BETA-1.96*curr_row$SE),2)},{round(exp(curr_row$BETA+1.96*curr_row$SE),2)})")
    }
  }
  
  # Label variants that are not present in NIAGADS HISP, to avoid having to recheck
  out_df %<>% mutate(NIA_P = ifelse(str_detect(ID,'1-143210077'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'1-15860133'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'11-55233749'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'17-75801070'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'18-33147978'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'19-48564935'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'20-36949029'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'1-143210080'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'12-79346706'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'4-154591668'),"ND",NIA_P)) %>%
    mutate(NIA_P = ifelse(str_detect(ID,'1-39179844'),"ND",NIA_P))
  
  # And other variants
  out_df %<>% mutate(AoU_AMR_P = ifelse(str_detect(ID,'17-75801070'),'ND',AoU_AMR_P))
  
  # hand off the output, clean it up a bit
  return(out_df %>% select(-BETA,-SE))
}
make_functional_table_aou = function(gw_hisp,gw_amr,gw_proj,gw_matchit) {
  # Make data frame
  out_df = data.frame(ID=gw_hisp$ID,CHR=gw_hisp$CHR,rsid=NA,Gene=NA,EnhancerSite=NA,EnhancerGene=NA,
                      CADD_Phred=NA,EpiActivMax=NA,EpiReprMax=NA,EpiTransMax=NA,Conserved=NA) %>%
    add_row(CADD_Phred = -1) %>%
    add_row(data.frame(ID=gw_amr$ID,CHR=gw_amr$CHR,rsid=NA,Gene=NA,EnhancerSite=NA,EnhancerGene=NA,
                       CADD_Phred=NA,EpiActivMax=NA,EpiReprMax=NA,EpiTransMax=NA,Conserved=NA)) %>%
    add_row(CADD_Phred = -1) %>%
    add_row(data.frame(ID=gw_proj$ID,CHR=gw_proj$CHR,rsid=NA,Gene=NA,EnhancerSite=NA,EnhancerGene=NA,
                       CADD_Phred=NA,EpiActivMax=NA,EpiReprMax=NA,EpiTransMax=NA,Conserved=NA)) %>%
    add_row(CADD_Phred = -1) %>%
    add_row(data.frame(ID=gw_matchit$ID,CHR=gw_matchit$CHR,rsid=NA,Gene=NA,EnhancerSite=NA,EnhancerGene=NA,
                       CADD_Phred=NA,EpiActivMax=NA,EpiReprMax=NA,EpiTransMax=NA,Conserved=NA))
  
  # Get favor_data
  favor_annot = vroom("../Data_Links/FAVOR/favor_hits_annot.txt",skip=1,delim=",",col_names=F)
  names(favor_annot) = names(vroom("../Data_Links/FAVOR/favor_hits_annot.txt",n_max=1,delim="\t"))
  print(glue("Nrow out_df {nrow(out_df %>% filter(!is.na(ID)))}. Overlap with FAVOR: {nrow(out_df %>% filter(ID %in% favor_annot$variant_vcf))}"))
  
  # Match annot to out_df
  super_enhancers = vroom("../00_AoU/dbSUPER_SuperEnhancers_hg19.tsv")
  for (row in 1:nrow(out_df)) {
    curr_row = out_df[row,]
    if (is.na(curr_row$ID)) next # deal with placeholder rows
    favor_row = favor_annot %>% filter(variant_vcf == curr_row$ID)
    out_df$Gene[[row]] = favor_row$genecode_comprehensive_info
    out_df$CADD_Phred[[row]] = round(favor_row$cadd_phred,2)
    out_df$EnhancerSite[[row]] = ifelse(!is.na(favor_row$cage_enhancer) |
                                            !is.na(favor_row$genehancer) |
                                            !is.na(favor_row$super_enhancer),
                                          "TRUE","FALSE")
    if (!is.na(favor_row$genehancer)) {
      lgenes = regmatches(favor_row$genehancer, gregexpr("(?<=connected_gene=)[^;]+",favor_row$genehancer, perl = TRUE))[[1]]
      out_df$EnhancerGene[[row]] = paste(lgenes,collapse=",")
    }
    if (!is.na(favor_row$super_enhancer)) {
      senhancers = str_split(favor_row$super_enhancer,",")[[1]]
      db_match = which(super_enhancers$se_id %in% senhancers)
      gene_matches = unique(super_enhancers$gene_symbol[db_match])
      out_df$EnhancerGene[[row]] = paste(gene_matches,collapse=",")
    }
    
    out_df$EpiActivMax[[row]] = round(max(c(favor_row$apc_epigenetics_active,favor_row$encodeh3k27ac_sum,
                                            favor_row$encodeh3k4me1_sum,favor_row$encodeh3k4me2_sum,
                                            favor_row$encodeh3k4me3_sum,favor_row$encodeh3k9ac_sum,
                                            favor_row$encodeh4k20me1_sum,favor_row$encodeh2afz_sum),na.rm=T),2)
    out_df$EpiReprMax[[row]] = round(max(c(favor_row$apc_epigenetics_repressed,favor_row$encodeh3k9me3_sum,
                                           favor_row$encodeh3k27me3_sum),na.rm=T),2)
    out_df$EpiTransMax[[row]] = round(max(c(favor_row$apc_epigenetics_transcription,
                                            favor_row$encodeh3k36me3_sum,
                                            favor_row$encodeh3k79me2_sum),na.rm=T),2)
    out_df$Conserved[[row]] = ifelse((!is.na(favor_row$priphylop) & favor_row$priphylop >= 0.3) |
                                       (!is.na(favor_row$priphcons) & favor_row$priphcons >= 0.3) |
                                       (!is.na(favor_row$mamphcons) & favor_row$mamphcons >= 0.3) |
                                       (!is.na(favor_row$verphcons) & favor_row$verphcons >= 0.3) |
                                       (!is.na(favor_row$gerp_n) & favor_row$gerp_n >= 10) |
                                       (!is.na(favor_row$gerp_s) & favor_row$gerp_s >= 10),
                                     "TRUE","FALSE")
  }
  # Clean gene name
  out_df %<>% mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
                     Gene = str_replace_all(Gene," ",""),
                     Gene = gsub("\\([^)]*\\)", "", Gene),
                     Gene = str_replace(Gene,'NONE,','')) %>%
    mutate(Gene = str_replace(Gene,'ENSG00000257894','Lnc-PAWR-1'))
  
  # Get rsids
  favor_rsids = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/FAVOR_Rsids.txt") %>%
    rename(ID = FAVOR_VCF)
  for (row in 1:nrow(out_df)) {
    match = favor_rsids %>% filter(ID == out_df$ID[[row]])
    if (nrow(match)>1) print(match)
    if (nrow(match)>0) out_df$rsid[[row]] = match$Rsid
  }
  
  # Remove variants that are less likely to have a functional impact
  unlikely_rows = which(is.na(out_df$EnhancerGene) & out_df$CADD_Phred < 10 & out_df$CADD_Phred != -1 &
                          out_df$EpiActivMax < 10 & out_df$EpiReprMax < 10 & 
                          out_df$EpiTransMax < 10 & out_df$Conserved == 'FALSE')
  out_df = out_df[-c(unlikely_rows),]
  
  return(out_df %>% filter(!is.na(CADD_Phred)) %>% mutate(CADD_Phred = ifelse(CADD_Phred == -1,NA,CADD_Phred)) %>%
           select(-EnhancerSite))
}
getLambda = function(dataset) {
  gc()
  if (dataset == "AoU_HISP_Hisp_PCs") {
    df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/HISP_GWAS_Hisp_PCs/aou_hisp_hisp_pcs_geno_1e-1_mac_20.txt")
    case = 1695 ; control = 54767
  }else if (dataset == "AoU_HISP_Multianc_PCs") {
    df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/HISP_GWAS_Multianc_PCs/aou_AD_any_grp_hisp_gwas.txt")
    case = 1695 ; control = 54767
  }else if (dataset == "AoU_NIA_PROJ_3SD") {
    df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AoU_NIAGADS_Proj_3SD/aou_Proj_3SD_summ_stats_AD_any.txt")
    case = 9218 ; control = 172569
  }else if (dataset == "AoU_AMR") {
    df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_amr_summ_stats_AD_any.txt")
    case = 1107 ; control = 43779
  }else if (dataset == "NIAGADS_HISP") {
    df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full.dn8.gz") %>%
      mutate(CHISQ = qchisq(1-p,1))
    case = 2624 ; control = 6060
  }
  lambda = median(df$CHISQ,na.rm=T)/qchisq(0.5,1)
  lambda1000 = 1 + (lambda - 1) * (1 / case + 1 / control) * 500
  out_vec = c(lambda,lambda1000)
  return(out_vec)
}
make_files_for_hwe = function() {
  # get CHRPOS of key variants for HWE testing. Not using IDS here as it is slower
  # The output also includes IDs, so I can do matching there.
  df_genpcs = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_NoAgeSex_No_Anc_Covar/aou_MatchIt_Gen_summ_stats_AD_any_p_1e-2.txt") %>%
    mutate(Pval = 10^(-LOG10P)) %>% filter(Pval <= 1e-5)
  df_genpcs_agesex = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_AgeSex_No_Anc_Covar/aou_MatchIt_summ_stats_AD_any_p_1e-1.txt") %>% 
    mutate(Pval = 10^(-LOG10P)) %>% filter(Pval <= 1e-5)
  df_nia_hisp = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.p_5e-2.txt") %>%
    filter(p <= 1e-5) %>% mutate(chr = as.numeric(chr))
  df_meta_aou_genpcs = vroom("/n/home09/jwillett/true_lab_storage/02_NIAGADS_HISP/meta/aou_nia_hisp_genmatch_metaanalysis_chrposrefalt_p_1e-2.TBL") %>% 
    filter(`P-value` <= 1e-5) %>% mutate(CHR = as.numeric(CHR))
  df_meta_aou_genpcs_phen = vroom("/n/home09/jwillett/true_lab_storage/02_NIAGADS_HISP/meta/aou_nia_hisp_genphenmatch_metaanalysis_chrposrefalt_p_1e-2.TBL") %>% 
    filter(`P-value` <= 1e-5) %>% mutate(CHR = as.numeric(CHR))
  
  bed_file = data.frame(CHR=df_genpcs$CHROM,POSS=df_genpcs$GENPOS,POSE=df_genpcs$GENPOS) %>%
    add_row(data.frame(CHR=df_genpcs_agesex$CHROM,POSS=df_genpcs_agesex$GENPOS,POSE=df_genpcs_agesex$GENPOS)) %>%
    add_row(data.frame(CHR=df_nia_hisp$chr,POSS=df_nia_hisp$pos,POSE=df_nia_hisp$pos)) %>%
    add_row(data.frame(CHR=df_meta_aou_genpcs$CHR,POSS=df_meta_aou_genpcs$POS,POSE=df_meta_aou_genpcs$POS)) %>%
    add_row(data.frame(CHR=df_meta_aou_genpcs_phen$CHR,POSS=df_meta_aou_genpcs_phen$POS,POSE=df_meta_aou_genpcs_phen$POS)) %>%
    distinct(.keep_all = T)
  vroom_write(bed_file,"working/all_study_hits.bed")
  print("Wrote BED file for AoU testing to: working/all_study_hits.bed")
}
get_meta_nominal_hits = function(dataset,p_cutoff=1e-5,window=500) {
  print("HWE data obtained similar to AoU UKB NIA paper, just using this study's IDs for P<1e-5")
  # FIRST, GET THE DATA FOR NOMINAL HITS
  if (dataset == "AoU_NIA_Gen_GW") {
    gw_df = vroom("meta/aou_nia_hisp_genmatch_metaanalysis_chrposrefalt_p_1e-2.TBL",show_col_types = F) %>%
      filter(`P-value` <= p_cutoff) %>% 
      arrange(CHR,POS) %>% mutate(nearestgene = NA) 
    hardy = vroom("working/anc_hwe_midp.txt") %>% 
      filter(MIDP_AFR > 1e-15,MIDP_AMR > 1e-15,MIDP_EAS > 1e-15, MIDP_EUR > 1e-15,
             MIDP_MID > 1e-15, MIDP_SAS > 1e-15) # remove AoU variants in HWE
    gw_df %<>% filter(ifelse(nchar(Allele1)==1 & nchar(Allele2)==1,
                             MarkerName %in% hardy$ID | MarkerName %in% hardy$IDrev,MarkerName %in% hardy$ID))
  }else if (dataset == "AoU_NIA_GenPhen_GW") {
    gw_df = vroom("meta/aou_nia_hisp_genphenmatch_metaanalysis_chrposrefalt_p_1e-2.TBL",show_col_types = F) %>%
      filter(`P-value` <= p_cutoff) %>% 
      arrange(CHR,POS) %>% mutate(nearestgene = NA) 
    hardy = vroom("working/anc_hwe_midp.txt") %>% 
      filter(MIDP_AFR > 1e-15,MIDP_AMR > 1e-15,MIDP_EAS > 1e-15, MIDP_EUR > 1e-15,
             MIDP_MID > 1e-15, MIDP_SAS > 1e-15) # remove AoU variants in HWE
    gw_df %<>% filter(ifelse(nchar(Allele1)==1 & nchar(Allele2)==1,
                             MarkerName %in% hardy$ID | MarkerName %in% hardy$IDrev,MarkerName %in% hardy$ID))
  }
  out_df = assignLociNumbers(gw_df %>% mutate(Locus=NA,NewOld = NA)) %>%
    filter(MaxFreq - MinFreq < 0.4)
  
  # SECOND, CLASSIFY HITS AS NOVEL VS OLD (USING CODE FROM AOU VS UKB META ANALYSIS)
  known_loci = vroom("/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/known_loci_chrpos.txt")
  for (row in 1:nrow(out_df)) {
    curr_b_loc = known_loci %>% filter(CHR == out_df$CHR[[row]],
                                       POS - window*1000 <= out_df$POS[[row]],
                                       POS + window*1000 >= out_df$POS[[row]]) 
    if (nrow(curr_b_loc) > 0) out_df$NewOld[[row]] = "Old"
    else out_df$NewOld[[row]] = "New"
  }
  
  # FOURTH, add CHRPOS column to make intersecting datasets easier
  out_df %<>% mutate(CHRPOS = glue("{CHR}-{POS}"))
  
  return(out_df %>% rename(Pval = `P-value`,ID=MarkerName))
}
get_favor_annot = function(meta_gen,meta_genphen) {
  out_df = data.frame(ID=NA,Gene=NA,CADD_Phred=NA,EnhancerRole=NA,EnhancerLinkedGene=NA,
                      EpiActivMax=NA,EpiReprMax=NA,EpiTransMax=NA,Conserved=NA,FAVORAnnot=NA,
                      AOUGENP=NA,AOUGENPHENP=NA,NIAP=NA,NewOld=NA)
  favor_annot = vroom("../Data_Links/FAVOR/favor_hits_annot.txt",skip=1,delim=",",col_names=F)
  names(favor_annot) = names(vroom("../Data_Links/FAVOR/favor_hits_annot.txt",n_max=1,delim="\t"))
  super_enhancers = vroom("../00_AoU/dbSUPER_SuperEnhancers_hg19.tsv")
  out_df %<>% add_row(data.frame(ID = meta_gen$ID,NewOld=meta_gen$NewOld)) %>% add_row() %>%
    add_row(data.frame(ID = meta_genphen$ID,NewOld=meta_genphen$NewOld))
  
  # GET PVAL ACROSS STUDIES
  aou_gen = vroom("working/meta_hits_gwas_genonly_intersect.txt",show_col_types = F) %>%
    mutate(ID=glue("{CHROM}-{GENPOS}-{ALLELE0}-{ALLELE1}"),IDrev=glue("{CHROM}-{GENPOS}-{ALLELE1}-{ALLELE0}"))
  aou_genphen = vroom("working/meta_hits_gwas_genphen_intersect.txt",show_col_types = F) %>%
    mutate(ID=glue("{CHROM}-{GENPOS}-{ALLELE0}-{ALLELE1}"),IDrev=glue("{CHROM}-{GENPOS}-{ALLELE1}-{ALLELE0}"))
  nia_hisp = vroom("working/meta_hits_gwas_nia_intersect.txt",show_col_types = F) %>%
    mutate(ID=glue("{chr}-{pos}-{otherallele}-{effectallele}")) %>%
    mutate(IDrev=glue("{chr}-{pos}-{effectallele}-{otherallele}"))
  for (row in 1:nrow(out_df)) {
    if (is.na(out_df$ID[[row]])) next
    curr_row=out_df[row,]
    aou_gen_match = aou_gen %>% filter(ID == curr_row$ID | IDrev == curr_row$ID)
    if (nrow(aou_gen_match)>0) out_df$AOUGENP[[row]] = aou_gen_match$Pval
    aou_genphen_match = aou_genphen %>% filter(ID == curr_row$ID | IDrev == curr_row$ID)
    if (nrow(aou_genphen_match)>0) out_df$AOUGENPHENP[[row]] = aou_genphen_match$Pval
    nia_hisp_match = nia_hisp %>% filter(ID == curr_row$ID | IDrev == curr_row$ID)
    if (nrow(nia_hisp_match)>0) out_df$NIAP[[row]] = nia_hisp_match$p
  }
  
  # GET FAVOR ANNOT
  for (row in 1:nrow(out_df)) {
    if (is.na(out_df$ID[[row]])) next
    favor_match = which(favor_annot$variant_vcf == out_df$ID[[row]] & !is.na(favor_annot$vid))
    favor_row = favor_annot[favor_match,]
    if (!is.na(favor_row$cadd_phred)) out_df$CADD_Phred[[row]] = favor_row$cadd_phred
    if (!is.na(favor_row$cage_enhancer) | !is.na(favor_row$genehancer) | !is.na(favor_row$super_enhancer)) {
      out_df$EnhancerRole[[row]] = "Yes"
      if (!is.na(favor_row$genehancer)) {
        # The score is a factor, but what is more important is the existence of a genehancer that provides
        # evidence that there is some relationship between the annotation and a gene's expression. The
        # scores are also not necessarily interpretable since they can diverge by so much.
        lgenes = regmatches(favor_row$genehancer, gregexpr("(?<=connected_gene=)[^;]+",favor_row$genehancer, perl = TRUE))[[1]]
        out_df$EnhancerLinkedGene[[row]] = paste(lgenes,collapse=",")
        # print("Gene linked to variant by GeneHancer")
      }
      if (!is.na(favor_row$super_enhancer)) {
        senhancers = str_split(favor_row$super_enhancer,",")[[1]]
        db_match = which(super_enhancers$se_id %in% senhancers)
        gene_matches = unique(super_enhancers$gene_symbol[db_match])
        out_df$EnhancerLinkedGene[[row]] = paste(gene_matches,collapse=",")
        # print("Gene linked to Variant by SuperEnhancer")
      }
    }
    
    out_df$EpiActivMax[[row]] = max(c(favor_row$apc_epigenetics_active,favor_row$encodeh3k27ac_sum,
                                      favor_row$encodeh3k4me1_sum,favor_row$encodeh3k4me2_sum,
                                      favor_row$encodeh3k4me3_sum,favor_row$encodeh3k9ac_sum,
                                      favor_row$encodeh4k20me1_sum,favor_row$encodeh2afz_sum),na.rm=T)
    out_df$EpiReprMax[[row]] = max(c(favor_row$apc_epigenetics_repressed,favor_row$encodeh3k9me3_sum,
                                     favor_row$encodeh3k27me3_sum),na.rm=T)
    out_df$EpiTransMax[[row]] = max(c(favor_row$apc_epigenetics_transcription,
                                      favor_row$encodeh3k36me3_sum,
                                      favor_row$encodeh3k79me2_sum),na.rm=T)
    
    # Full Favor annotation
    favor_sig = ""
    if (!is.na(favor_row$cage_enhancer)) favor_sig %<>% paste("CAGE Enhancer;")
    if (!is.na(favor_row$genehancer)) favor_sig %<>% paste(glue("Genehancer;"))
    if (!is.na(favor_row$super_enhancer)) favor_sig %<>% paste("SuperEnhancer;")
    if (!is.na(favor_row$apc_conservation) & favor_row$apc_conservation >= 10) favor_sig %<>% paste(glue("aPC-Conservation {favor_row$apc_conservation};"))
    if (!is.na(favor_row$apc_epigenetics_active) & favor_row$apc_epigenetics_active >= 10) favor_sig %<>% paste(glue("aPC-Epigenetics-Active {favor_row$apc_epigenetics_active};"))
    if (!is.na(favor_row$apc_epigenetics_transcription) & favor_row$apc_epigenetics_transcription >= 10) favor_sig %<>% paste(glue("aPC-Epigenetics-Transcription {favor_row$apc_epigenetics_transcription};"))
    if (!is.na(favor_row$apc_epigenetics_repressed) & favor_row$apc_epigenetics_repressed >= 10) favor_sig %<>% paste(glue("aPC-Epigenetics-Repressed {favor_row$apc_epigenetics_repressed};"))
    if (!is.na(favor_row$apc_local_nucleotide_diversity) & favor_row$apc_local_nucleotide_diversity >= 10) favor_sig %<>% paste(glue("aPC-Local-Nucleotide-Diversity {favor_row$apc_local_nucleotide_diversity};"))
    if (!is.na(favor_row$apc_transcription_factor) & favor_row$apc_transcription_factor >= 10) favor_sig %<>% paste(glue("aPC-Transcription-Factor {favor_row$apc_transcription_factor};"))
    if (!is.na(favor_row$apc_mappability) & favor_row$apc_mappability >= 10) favor_sig %<>% paste(glue("aPC-Mappability {favor_row$apc_mappability};"))
    if (!is.na(favor_row$apc_mutation_density) & favor_row$apc_mutation_density >= 10) favor_sig %<>% paste(glue("aPC-Mutation-Density {favor_row$apc_mutation_density};"))
    if (!is.na(favor_row$cadd_phred) & favor_row$cadd_phred >= 10) favor_sig %<>% paste(glue("CADD Phred {favor_row$cadd_phred};"))
    if (!is.na(favor_row$encode_dnase_sum) & favor_row$encode_dnase_sum >= 0.44) favor_sig %<>% paste(glue("Active DNase {favor_row$encode_dnase_sum};"))
    if (!is.na(favor_row$encodeh2afz_sum) & favor_row$encodeh2afz_sum >= 3.28) favor_sig %<>% paste(glue("Active H2AFZ {favor_row$encodeh2afz_sum};"))
    if (!is.na(favor_row$encodeh3k4me1_sum) & favor_row$encodeh3k4me1_sum >= 4.5) favor_sig %<>% paste(glue("Active H3K4me1 {favor_row$encodeh3k4me1_sum};"))
    if (!is.na(favor_row$encodeh3k4me2_sum) & favor_row$encodeh3k4me2_sum >= 3.5) favor_sig %<>% paste(glue("Active H3K4me2 {favor_row$encodeh3k4me2_sum};"))
    if (!is.na(favor_row$encodeh3k4me3_sum) & favor_row$encodeh3k4me3_sum >= 3.7) favor_sig %<>% paste(glue("Active H3K4me3 {favor_row$encodeh3k4me3_sum};"))
    if (!is.na(favor_row$encodeh3k9ac_sum) & favor_row$encodeh3k9ac_sum >= 3.1) favor_sig %<>% paste(glue("Active H3K9ac {favor_row$encodeh3k9ac_sum};"))
    if (!is.na(favor_row$encodeh3k9me3_sum) & favor_row$encodeh3k9me3_sum >= 3.7) favor_sig %<>% paste(glue("Repressed H3K9me3 {favor_row$encodeh3k9me3_sum};"))
    if (!is.na(favor_row$encodeh3k27ac_sum) & favor_row$encodeh3k27ac_sum >= 4.5) favor_sig %<>% paste(glue("Active H3K27ac {favor_row$encodeh3k27ac_sum};"))
    if (!is.na(favor_row$encodeh3k27me3_sum) & favor_row$encodeh3k27me3_sum >= 4.5) favor_sig %<>% paste(glue("Repressed H3K27me3 {favor_row$encodeh3k27me3_sum};"))
    if (!is.na(favor_row$encodeh3k36me3_sum) & favor_row$encodeh3k36me3_sum >= 3.9) favor_sig %<>% paste(glue("Transcription H3K36me3 {favor_row$encodeh3k36me3_sum};"))
    if (!is.na(favor_row$encodeh3k79me2_sum) & favor_row$encodeh3k79me2_sum >= 3.5) favor_sig %<>% paste(glue("Transcription H3k79me2 {favor_row$encodeh3k79me2_sum};"))
    if (!is.na(favor_row$encodeh4k20me1_sum) & favor_row$encodeh4k20me1_sum >= 3.7) favor_sig %<>% paste(glue("Active H4k20me1 {favor_row$encodeh4k20me1_sum};"))
    if (!is.na(favor_row$encodetotal_rna_sum) & favor_row$encodetotal_rna_sum >= 0.1) favor_sig %<>% paste(glue("Transcription totalRNA {favor_row$encodetotal_rna_sum};"))
    if (!is.na(favor_row$priphylop) & favor_row$priphylop >= 0.3) favor_sig %<>% paste(glue("Conservation priPhyloP {favor_row$priphylop};"))
    if (!is.na(favor_row$priphcons) & favor_row$priphcons >= 0.3) favor_sig %<>% paste(glue("Conservation priPhCons {favor_row$priphcons};"))
    if (!is.na(favor_row$mamphcons) & favor_row$mamphcons >= 0.3) favor_sig %<>% paste(glue("Conservation mamPhCons {favor_row$mamphcons};"))
    if (!is.na(favor_row$verphcons) & favor_row$verphcons >= 0.3) favor_sig %<>% paste(glue("Conservation verPhCons {favor_row$verphcons};"))
    if (!is.na(favor_row$gerp_n) & favor_row$gerp_n >= 10) favor_sig %<>% paste(glue("Conservation GerpN {favor_row$gerp_n};"))
    if (!is.na(favor_row$gerp_s) & favor_row$gerp_s >= 10) favor_sig %<>% paste(glue("Conservation GerpS {favor_row$gerp_s};"))
    
    out_df$Gene[[row]] = favor_row$genecode_comprehensive_info
    out_df$FAVORAnnot[[row]] = favor_sig
    if (str_detect(favor_sig,"Conservation")) out_df$Conserved[[row]] = "Yes"
    else out_df$Conserved[[row]] = "No"
  }
  return(out_df)
}
make_favor_functional_table = function(df) {
  out_df = df
  rows_to_cut = numeric()
  
  for (row in 1:nrow(df)) {
    if (is.na(df$ID[[row]]) | row %in% rows_to_cut) next
    curr_row = df[row,]
    
    # FUNCTIONAL ROLE?
    if (curr_row$CADD_Phred < 10 & is.na(curr_row$EnhancerRole) &
        curr_row$EpiActivMax < 10 & curr_row$EpiReprMax < 10 &
        curr_row$EpiTransMax < 10 & curr_row$Conserved == "No")
      rows_to_cut %<>% append(row)
    
    # NOMINAL IN BOTH DATASETS?
    if (curr_row$AOUGENP > 0.05 & curr_row$AOUGENPHENP > 0.05) rows_to_cut %<>% append(row)
    else if (curr_row$NIAP > 0.05) rows_to_cut %<>% append(row)
  }
  return(out_df[-c(rows_to_cut),])
}
get_gene_expression = function(df,p_cutoff=0.01) {
  out_df = df %>% select(ID,CHR,Rsid,Gene,MAF) %>%
    mutate(Cell=NA)
  sc_results_general = vroom("../00_AoU/scrna_logfold_stats_by_cell_group_by_pathology_and_symptoms_1234.txt") %>%
    filter(PValAdj <= p_cutoff) %>% arrange(Gene,desc(AvgLog2FC)) %>% group_by(Gene) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"second_vs_first"),"2v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"third_vs_first"),"3v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_first"),"4v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"third_vs_second"),"3v2",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_second"),"4v2",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_third"),"4v3",Comparison)) 
  
  for (row in 1:nrow(out_df)) {
    genes = unique(str_split(out_df$Gene[[row]],",")[[1]])
    clean_str = ""

    for (gene in genes) { # to deal with more than one gene
      gene_edit = str_replace_all(gene,"\\(dist.*?\\)", "")
      gene_edit = str_replace_all(gene_edit," ","")
      gene_edit = gsub("\\([^)]*\\)", "", gene_edit)
      cell_match = sc_results_general %>% filter(Gene == gene_edit)
      if (nrow(cell_match)<1) next
      
      cells = c("Exc ","Inh ","Oli ","OPC ","Ast ","Mic ")
      comps = c("2v1","3v1","4v1","3v2","4v2","4v3")
      
      for (cell in cells) {
        for (comp in comps) {
          cell_match_rows = cell_match %>% filter(str_detect(CellPop,cell),str_detect(Comparison,comp)) %>%
            mutate(Sign=ifelse(AvgLog2FC >= 0,"+","-"))
          count = nrow(cell_match_rows)
          if (count > 0) {
            sign = paste0(cell_match_rows$Sign,collapse = "")
            clean_str = glue("{clean_str} \\({gene_edit}\\) {cell} {comp} (x{count}) {sign}")
          }
        }
      }
    }
    out_df$Cell[[row]] = clean_str
    if (length(genes) == 1) {
      out_df$Cell[[row]] = str_replace_all(out_df$Cell[[row]],glue("\\({genes[[1]]}\\) "),"")
    }
  }
  return(out_df)
}
get_pathspecific_gene_expression = function(df,p_cutoff = 0.01) {
  out_df = data.frame(ID=NA,CHR=NA,Rsid=NA,Gene=NA,MAF=NA,Path=NA,Cell=NA)
  sc_results = vroom("../00_AoU/scrna_logfold_stats_by_cell.txt") %>%
    filter(PValAdj <= p_cutoff) %>% arrange(Gene,desc(AvgLog2FC)) %>% 
    mutate(Comparison = str_replace(Comparison,".csv","")) %>% group_by(Gene) %>%
    mutate(CellPop = ifelse(str_detect(Comparison,"second_vs_first"),glue("{CellPop} 2v1"),CellPop)) %>%
    mutate(CellPop = ifelse(str_detect(Comparison,"third_vs_second"),glue("{CellPop} 3v2"),CellPop))
  
  for (row in 1:nrow(df)) {
    for (path in unique(sc_results$Outcome)) {
      genes = unique(str_split(df$Gene[[row]],",")[[1]])
      clean_str = ""
      
      for (gene in genes) { # to deal with more than one gene
        gene_edit = str_replace_all(gene,"\\(dist.*?\\)", "")
        gene_edit = str_replace_all(gene_edit," ","")
        gene_edit = gsub("\\([^)]*\\)", "", gene_edit)
        cell_match = sc_results %>% filter(Gene == gene_edit,Outcome == path)
        if (nrow(cell_match)<1) next
        
        cells = c("Exc ","Inh ","Oli ","OPC ","Ast ","Mic ")
        comps = c("2v1","3v1","3v2")
        
        for (cell in cells) {
          for (comp in comps) {
            cell_match_rows = cell_match %>% filter(str_detect(CellPop,cell),str_detect(CellPop,comp)) %>%
              mutate(Sign=ifelse(AvgLog2FC >= 0,"+","-"))
            count = nrow(cell_match_rows)
            if (count > 0) {
              sign = paste0(cell_match_rows$Sign,collapse = "")
              clean_str = glue("{clean_str} ({gene_edit}) {cell} {comp} (x{count}) {sign}")
            }
          }
        }
      }
      out_df %<>% add_row(data.frame(ID=df$ID[[row]],CHR=df$CHR[[row]],Rsid=df$Rsid[[row]],
                                     Gene=df$Gene[[row]],MAF=df$MAF[[row]],
                                     Path=path,Cell=clean_str))
    }
  }
  return(out_df)
}
