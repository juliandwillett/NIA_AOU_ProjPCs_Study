### Check stats for GW-sig hits
ancestries=(afr amr eas eur mid sas) ;\
for anc in "${ancestries[@]}"; do \
  for ((chr=1;chr<=22;chr++)); do \
    # AoU HISP
    ./plink2 --pfile pgen_qc/chr${chr}_geno_mac --extract aou_hisp_for_hardy.txt \
      --keep ids/hispanic_individuals.txt ${anc}_ids.txt --hardy 'midp' --out variant_qc/aou_hisp_chr${chr}_anc_${anc}
  
    # AoU AMR
    ./plink2 --pfile pgen_qc/chr${chr}_geno_mac --extract aou_amr_for_hardy.txt \
      --keep ids/amr_ids.txt --hardy 'midp' --out variant_qc/aou_amr_chr${chr}_anc_amr
  
    # AoU NIA PROJ 3SD
    ./plink2 --pfile pgen_qc/chr${chr}_geno_mac --extract aou_nia_proj_3sd_for_hardy.txt \
      --keep piezo2_projection/aou_within_3sd_projected_niagads_ids.txt ${anc}_ids.txt \
      --hardy 'midp' --out variant_qc/aou_nia_proj_3sd_chr${chr}_anc_${anc}

    # AoU NIA PROJ MatchIt
    ./plink2 --pfile pgen_qc/chr${chr}_geno_mac --extract aou_nia_matchit_for_hardy.txt \
      --keep ids/ids_matchit.txt ${anc}_ids.txt \
      --hardy 'midp' --out variant_qc/aou_nia_matchit_chr${chr}_anc_${anc}
  done
done

### Merge the output 
# HWE
for anc in "${ancestries[@]}"; do \
  echo $anc ;\
  head -1 variant_qc/aou_hisp_chr1_anc_${anc}.hardy > variant_qc/aou_hisp_all_hardy_anc_${anc}.txt ;\
  head -1 variant_qc/aou_amr_chr1_anc_amr.hardy > variant_qc/aou_amr_all_hardy_anc_amr.txt ;\
  head -1 variant_qc/aou_nia_proj_3sd_chr1_anc_${anc}.hardy > variant_qc/aou_nia_proj_3sd_all_hardy_anc_${anc}.txt ;\
  head -1 variant_qc/aou_nia_matchit_chr1_anc_${anc}.hardy > variant_qc/aou_nia_matchit_all_hardy_anc_${anc}.txt ;\
  for ((chr=1;chr<=22;chr++)); do \
    tail -n +2 variant_qc/aou_hisp_chr${chr}_anc_${anc}.hardy >> variant_qc/aou_hisp_all_hardy_anc_${anc}.txt ;\
    tail -n +2 variant_qc/aou_amr_chr${chr}_anc_amr.hardy >> variant_qc/aou_amr_all_hardy_anc_amr.txt ;\
    tail -n +2 variant_qc/aou_nia_proj_3sd_chr${chr}_anc_${anc}.hardy >> variant_qc/aou_nia_proj_3sd_all_hardy_anc_${anc}.txt ;\
    tail -n +2 variant_qc/aou_nia_matchit_chr${chr}_anc_${anc}.hardy >> variant_qc/aou_nia_matchit_all_hardy_anc_${anc}.txt ;\
  done \
done

#### Then merge into a dataframe for all ancestries in R
{R}
library(vroom)
library(tidyverse)
library(magrittr)
library(glue)

datasets = c('hisp','amr','nia_proj_3sd','nia_matchit')
ancestries = c('afr','amr','eas','eur','mid','sas')
for (dataset in datasets) {
    df = data.frame(ID = read.table(glue("aou_{dataset}_for_hardy.txt"))$V1)
    for (anc in ancestries) {
        if (dataset == 'amr' & anc != 'amr') next
        hwe = vroom(glue("variant_qc/aou_{dataset}_all_hardy_anc_{anc}.txt"),show_col_types = F) %>% 
            select(ID,MIDP)
        df = merge(df,hwe,by='ID') 
        names(df)[[length(df)]] = glue("MIDP_{anc}")
    }
    vroom_write(df,glue("variant_qc/aou_{dataset}_all_hardy_anc_all.txt"))
}
{/R}
