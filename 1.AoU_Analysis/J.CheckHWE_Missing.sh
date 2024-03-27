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
