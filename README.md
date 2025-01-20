# Projected Principal Components Analysis

This is the workspace for our manuscript, "Matching Heterogeneous Cohorts by Projected Principal Components Reveals Two Novel Alzheimer's Disease-Associated Genes in the Hispanic Population."

To replicate our analyses, follow the flow as documented in the paper. So:
1. Organize data for first cohort (make bed files), running GWAS. 
2. Run the AoU_Analysis folder code to produce projected PCs onto the second cohort.
3. Run regenie on your second cohort using the projected PCs as covariates. 
4. Meta-analyze your second cohort with the first cohort using their respective GWAS files. 

FAQ:
Q: How did you filter by HWE?
A: We used a similar process of filtering by HWE by single ancestry, which we document in another paper (where the code is readily available). The paper has been accepted to Alzheimer's and Dementia. https://www.medrxiv.org/content/10.1101/2024.09.11.24313439v1
