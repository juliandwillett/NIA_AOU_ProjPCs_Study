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
  awk 'NR==1 {$5 = "ALLELE1"; $6 = "ALLELE0"; $7 = "A1FREQ"; $8 = "BETA"; $9 = "SE"; $10 = "Pval"; print $0} NR>1 {print $0}' \
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
PROCESS /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AoU_NIA_HISP_MatchIt_No_Anc_Covar/aou_MatchIt_summ_stats_AD_any_pval.txt
PROCESS /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfHISP8467_a_pc_5JPCs_ss.Affection.Status.glm.logistic.full_for_metal.txt

ANALYZE HETEROGENEITY

# For running that file:
salloc -p test --mem 80000 -t 0-02:00 -n 4
/n/home13/dprokopenko/bin/metal METAL_script_AOU_NIAGADS.txt

### Then rename output for downstream work.
mv 
