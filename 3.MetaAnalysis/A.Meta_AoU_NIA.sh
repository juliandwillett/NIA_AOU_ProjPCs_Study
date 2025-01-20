### format data, then run METAL.
aou_nia_genphen_match="/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_No_Anc_Covar/aou_MatchIt_summ_stats_AD_any" #.txt
nia_hisp="/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full" #.dn8.gz

### Match column names, as relevant for METAL
#MARKER   ID
#ALLELE   ALLELE1 ALLELE0
#FREQ     A1FREQ
#EFFECT   BETA
#STDERR   SE
#PVAL     Pval

awk 'NR==1 {$14 = "Pval"; print $0} NR>1 {$14 = 10^(-1 * $12); print $0}' ${aou_nia_genphen_match}.txt > ${aou_nia_genphen_match}_pval.txt
zcat ${nia_hisp}.dn8.gz |\
  awk 'NR==1 {$3 = "ID"; $5 = "ALLELE1"; $6 = "ALLELE0"; $7 = "A1FREQ"; $8 = "BETA"; $9 = "SE"; $10 = "Pval"; print $0} \
    NR>1 {gsub(":", "-", $3); print $0}' \
  > ${nia_hisp}_for_metal.txt

### Then run the meta analysis
# METAL_script_AOU_NIAGADS.txt
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON

MARKER   ID
ALLELE   ALLELE1 ALLELE0
FREQ     A1FREQ
EFFECT   BETA
STDERR   SE
PVAL     Pval
PROCESS /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_PCs_AgeSex_No_Anc_Covar/aou_MatchIt_summ_stats_AD_any_pval.txt
PROCESS /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full_for_metal.txt

ANALYZE HETEROGENEITY

# For running that file:
/n/home13/dprokopenko/bin/metal METAL_script_AOU_NIAGADS_genonly.txt
mv METAANALYSIS1.TBL meta/aou_nia_hisp_genmatch_metaanalysis.TBL
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0 "\tCHR\tPOS\tIDrev"} NR>1 {split($1, values, "-"); $16 = values[1]; \
  $17 = values[2]; $18 = $16 "-" $17 "-" toupper($3) "-" toupper($2); $20 = $16 "-" $17; print $0}' \
  meta/aou_nia_hisp_genmatch_metaanalysis.TBL > meta/aou_nia_hisp_genmatch_metaanalysis_chrposrefalt.TBL
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0} NR>1 && $10 <= 5e-8 && $16 != 23 {print $0}' \
  meta/aou_nia_hisp_genmatch_metaanalysis_chrposrefalt.TBL > \
  meta/aou_nia_hisp_genmatch_metaanalysis_chrposrefalt_gwsig.TBL

# FOR GEN PCs AND AGE AND SEX ANALYSIS
/n/home13/dprokopenko/bin/metal METAL_script_AOU_NIAGADS_genphen.txt
mv METAANALYSIS1.TBL meta/aou_nia_hisp_genphenmatch_metaanalysis.TBL
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0 "\tCHR\tPOS\tIDrev"} NR>1 {split($1, values, "-"); $16 = values[1]; \
  $17 = values[2]; $18 = $16 "-" $17 "-" toupper($3) "-" toupper($2); $20 = $16 "-" $17; print $0}' \
  meta/aou_nia_hisp_genphenmatch_metaanalysis.TBL > meta/aou_nia_hisp_genphenmatch_metaanalysis_chrposrefalt.TBL
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0} NR>1 && $10 <= 5e-8 && $16 != 23 {print $0}' \
  meta/aou_nia_hisp_genphenmatch_metaanalysis_chrposrefalt.TBL > \
  meta/aou_nia_hisp_genphenmatch_metaanalysis_chrposrefalt_gwsig.TBL
