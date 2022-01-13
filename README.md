# NicheSVM

NicheSVM is a tool for analyzing PICs (Physically Interacting Cells) to detect niche-specific expression.
Step1 - Calculating z-values
Step2 - Finding cell type DEGs and drawing heatmaps
Step3 - Deconvolution *Note that you should choose the clusters with unique DEGs (see the figure heatmap_clusterOnly~.pdf)
Step4 - Finding neighbor-specific DEGs

## Example

1. runNicheSVM_mouseEmbryoE7_5.m - PIC-seq data from mouse Embryo at E7.5
2. runNicheSVM_mouseEmbryoE8_5.m - PIC-seq data from mouse Embryo at E8.5
3. runNicheSVM_mouseEmbryoE9_5.m - PIC-seq data from mouse Embryo at E9.5
4. runNicheSVM_mouseHippo.m - Slide-seq data (Puck_200115_08)
5. runNicheSVM_liverPairedNormalized.m - pcRNA-seq data from mouse liver
