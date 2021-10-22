%%%%%%%%%%%%% Load MATLAB file %%%%%%%%%%%%%
load 'mouseEmbryoPIC.mat'
% 1. 'cell_id' - 33034 single cells + PICs
% 2. 'gene_name' - 2000 genes
% 3. 'log_data' - 2000 by 33034
% 4. 'clustering6color' - 33034 cells
% 5. 'clustering6name' - 33034 cells
% 6. 'clustering6name_unique' - 6 cell types
% 7. 'clustering8color' - 33034 cells
% 8. 'clustering8name' - 33034 cells
% 9. 'clustering8name_unique' - 8 cell types
% 10. 'clustering13color' - 33034 cells
% 11. 'clustering13name' - 33034 cells
% 12. 'clustering13name_unique' - 13 cell types
% 13. 'singletsIndex' - 33034 Boolean memberships for single cells
% 13. 'doubletsIndex' - 33034 Boolean memberships for PICs

%%%%%%%%%%%%% Calculating z-value %%%%%%%%%%%%%
log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
log_data_doublets=log_data_total(:,doubletsIndex);
log_data_doublets_zvalue=(log_data_doublets-repmat(mean(log_data_doublets,2),1,size(log_data_doublets,2)))./repmat(std(log_data_doublets')',1,size(log_data_doublets,2));

log_data_zvalue(isnan(log_data_zvalue))=0;
log_data_doublets_zvalue(isnan(log_data_doublets_zvalue))=0;

%%%%%%%%%%%%% ---- Pipeline for clustering13 ---- %%%%%%%%%%%%%
folderName='mouseEmbryo';
%%%%%%%%%%%%% 1) DEG by clustering 13 %%%%%%%%%%%%%
clusterSize=max(clustering13color);
[pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data,clustering13color);
save([folderName,'/pvalue_fdr_logRatio_zvalue.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')

%%%%%%%%%%%%% 2) PIC SVM classification %%%%%%%%%%%%%
clusterSelect=[1,2,3,5,6,7,8,9,10,11,12,13];
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
pCutoff=0.01;lrCutoff=0.4;
seedNumber=1;randSize=10000;
DEGnumber=5;
[bestMatch,artificialDoubletsCombiUnique]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering13color,clusterSelect,clustering13name_unique,log_data_zvalue,log_data_doublets_zvalue,DEGnumber);
save([folderName,'/SVM_bestMatch.mat'],'bestMatch','artificialDoubletsCombiUnique')

%%% Draw heatmap
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
clusterSelect=[1,2,3,5,6,7,8,9,10,11,12,13];
outputFile=[folderName,'/heatmap_PICSVM_top5DEG.pdf'];
pCutoff=0.01;lrCutoff=0.4;
DEGnumber=5;
drawHeatmap_PICSVM(outputFile,bestMatch,artificialDoubletsCombiUnique,clusterSelect,clustering13name_unique,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,gene_name,DEGnumber);

%%%%%%%%%%%%% 3) DEG between PICs and artificial doublets %%%%%%%%%%%%%
clusterSelect=[1,2,3,5,6,7,8,9,10,11,12,13];
seedNumber=1;randSize=10000;
minCell=5;
[pvalue_total,fdr_total,logRatio_total]=DEG_PIC_vs_AD(bestMatch,seedNumber,randSize,clustering13color,clusterSelect,log_data_zvalue,clustering13name_unique,log_data_doublets_zvalue,minCell);
save([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'],'pvalue_total','fdr_total','logRatio_total')

outputFolder=[folderName,'/heatmapPIC_vs_AD'];
load([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'])
pvalue_totalPIC_AD=pvalue_total;
logRatio_totalPIC_AD=logRatio_total;
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
pCutoff=0.01;lrCutoff=0.4;
DEGnumber=5;DEGnumberPIC_AD=10;
drawHeatmap_PIC_vs_AD(outputFolder,bestMatch,artificialDoubletsCombiUnique,clustering13color,clusterSelect,clustering13name_unique,log_data_zvalue,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name,DEGnumber,DEGnumberPIC_AD);

