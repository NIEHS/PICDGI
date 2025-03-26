# PICDGI: Predict Immunosuppressive Cancer Driver Gene using gene-gene Interaction feature

$~~$

**Traditional frequency-based methods identify driver mutations by evaluating their recurrence across tumors, operating under the assumption that mutations occurring more frequently are functionally significant. However, these approaches often overlook rare driver mutations, particularly those with low mutation frequencies.**

**To overcome this limitation, we present PICDGI, a computational framework designed to predict cancer driver genes (CDGs) by leveraging single-cell RNA sequencing (scRNA-seq) data to model gene-gene interactions, thereby capturing cellular heterogeneity within tumors.**

**Unlike conventional methods, PICDGI utilizes variational Bayesian inference to model gene expression dynamics and quantify the influence of gene-gene interactions on tumor progression. It identifies CDGs by detecting their immunosuppressive effects on immune cells through interaction-driven mechanisms during cancer development.**

**Applied to LUAD scRNA-seq data, PICDGI identified top-ranked CDGs, with 62% overlapping known OGs and TSGs, while 38% are novel candidates. Compared to Moran’s I test in Monocle3, PICDGI-identified genes showed higher expression, reinforcing their role in tumor growth. By capturing gene expression dynamics and interactions, PICDGI offers insights into potential therapeutic targets for personalized cancer treatment.**

![](Figure/EDO/picdgi/final/github.png)
