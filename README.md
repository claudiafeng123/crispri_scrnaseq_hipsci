# crispri_scrnaseq_hipsci
Analyzing data from two large-scale CRISRPi screens in iPSCs with single-cell RNA-seq read-out.

# Pipeline

1. Run CellRanger
2. Assign cells to donors (Vireo)
3. Quality control
4. Guide and donor assignment.
5. a) Merge per-inlet counts into a single object. b) UMAPs. c) Variance decomposition of expression. d) Control vs. unassigned cell comparison.
6. Compute expression fold-change compared to control cells.  
7. Experiment summary.
8. Reproducibility across guides, timepoints, experiments. Power estimation.
9. Target gene down-regulation per-line.
10. Trans effect analysis.
11. Target-target correlation analysis.
12. Expressed gene-expressed gene correlation analysis.
13. Pluripotency networks
14. Variation due to CRISPRi, guide, cell line and donor.
15. Heritability estimates.
16. Variance explained due to known cis eQTLs.
