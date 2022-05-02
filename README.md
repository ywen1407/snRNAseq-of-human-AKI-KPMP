# snRNAseq-of-human-AKI-KPMP

This is the code for manuscript titled "Analysis of the Human Kidney Transcriptome and Plasma Proteome Identifies Novel Biomarkers of Proximal Tubule Maladaptation to Injury". 

Author lists:

Authors:  Yumeng Wen1, Emily Su2, Leyuan Xu3, Steven Menez1, Dennis Moledina3, Paul M. Palevsky4, Lloyd Cantley3, Patrick Cahan2, Chirag R. Parikh1*, for the Kidney Precision Medicine Project (KPMP) and Translational Investigation of Biomarker Endpoint of Acute Kidney Injury (TRIBE-AKI) Consortia

Affiliations:
1. Department of Medicine/Division of Nephrology, Johns Hopkins University School of Medicine, Baltimore, MD, USA.
2. Department of Biomedical Engineering, Johns Hopkins University School of Medicine, Baltimore, MD, USA.
3. Department of Medicine/Section of Nephrology, Yale School of Medicine, New Haven, CT, USA.
4. Renal-Electrolyte Division, University of Pittsburgh School of Medicine, Pittsburgh, PA, USA.
*Corresponding author. Email:  chirag.parikh@jhmi.edu; Telephone: 410-955-5268

All dataset can be accessed via the publicly available KPMP respository at: atlas.kpmp.org/repository

The general work flow is:
1. decontX to remove contamination
2. snRNA-seq dataset preprocess, remove potential doublets, re-preprocess, and visualize cannonical marker gene expression among clusters of kidney cells. 
3. differential gene expression and FGSEA
4. human vs. mouse comparision (2 annotations): global and differentiating proximal tubular cells at differnet health states
- pysinglecellnet
- support vector machine
- SCMAP
5. gene regulatory network analysis with Epoch and pySCENIC

(to be continued) 
