### Check stats for GW-sig hits
for ((chr=1;chr<=22;chr++)); do 
  # AoU HISP
  ./plink2 --pfile pgen_qc/chr${chr}_geno_mac --extract variant_qc/gw_hits_aou_hisp.txt \
    --keep ids/hispanic_individuals.txt --hardy 'midp' --missing --out variant_qc/aou_hisp_chr${chr}

  # AoU AMR
  ./plink2 --pfile pgen_qc/chr${chr}_geno_mac --extract variant_qc/gw_hits_aou_amr.txt \
    --keep ids/amr_ids.txt --hardy 'midp' --missing --out variant_qc/aou_amr_chr${chr}

  # AoU NIA PROJ 3SD
  ./plink2 --pfile pgen_qc/chr${chr}_geno_mac --extract variant_qc/gw_hits_aou_nia_proj.txt \
    --keep piezo2_projection/aou_within_3sd_projected_niagads_ids.txt \
    --hardy 'midp' --missing --out variant_qc/aou_nia_proj_chr${chr}
done

### Merge the output
# VMISS
head -1 variant_qc/aou_hisp_chr1.vmiss > variant_qc/aou_hisp_all_vmiss.txt ;\
head -1 variant_qc/aou_amr_chr1.vmiss > variant_qc/aou_amr_all_vmiss.txt ;\
head -1 variant_qc/aou_nia_proj_chr1.vmiss > variant_qc/aou_nia_proj_all_vmiss.txt ;\
for ((chr=1;chr<=22;chr++)); do \
  tail -n +2 variant_qc/aou_hisp_chr${chr}.vmiss >> variant_qc/aou_hisp_all_vmiss.txt ;\
  tail -n +2 variant_qc/aou_amr_chr${chr}.vmiss >> variant_qc/aou_amr_all_vmiss.txt ;\
  tail -n +2 variant_qc/aou_nia_proj_chr${chr}.vmiss >> variant_qc/aou_nia_proj_all_vmiss.txt ;\
done

# SMISS
head -1 variant_qc/aou_hisp_chr1.smiss > variant_qc/aou_hisp_all_smiss.txt ;\
head -1 variant_qc/aou_amr_chr1.smiss > variant_qc/aou_amr_all_smiss.txt ;\
head -1 variant_qc/aou_nia_proj_chr1.smiss > variant_qc/aou_nia_proj_all_smiss.txt ;\
for ((chr=1;chr<=22;chr++)); do \
  tail -n +2 variant_qc/aou_hisp_chr${chr}.smiss >> variant_qc/aou_hisp_all_smiss.txt ;\
  tail -n +2 variant_qc/aou_amr_chr${chr}.smiss >> variant_qc/aou_amr_all_smiss.txt ;\
  tail -n +2 variant_qc/aou_nia_proj_chr${chr}.smiss >> variant_qc/aou_nia_proj_all_smiss.txt ;\
done

# HWE
head -1 variant_qc/aou_hisp_chr1.hardy > variant_qc/aou_hisp_all_hardy.txt ;\
head -1 variant_qc/aou_amr_chr1.hardy > variant_qc/aou_amr_all_hardy.txt ;\
head -1 variant_qc/aou_nia_proj_chr1.hardy > variant_qc/aou_nia_proj_all_hardy.txt ;\
for ((chr=1;chr<=22;chr++)); do \
  tail -n +2 variant_qc/aou_hisp_chr${chr}.hardy >> variant_qc/aou_hisp_all_hardy.txt ;\
  tail -n +2 variant_qc/aou_amr_chr${chr}.hardy >> variant_qc/aou_amr_all_hardy.txt ;\
  tail -n +2 variant_qc/aou_nia_proj_chr${chr}.hardy >> variant_qc/aou_nia_proj_all_hardy.txt ;\
done
