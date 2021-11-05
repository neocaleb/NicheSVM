%%%%%%%%%%%%% Load MATLAB file %%%%%%%%%%%%%
load 'mouseEmbryoPIC.mat'
% 1. 'cell_id' - 6288 single cells
% 2. 'cell_id_total' - 6288 single cells + 367 PICs
% 3. 'gene_name' - 1608 genes
% 4. 'log_data' - 1608 by 6288 single cells
% 5. 'log_data_total' - 1608 by 6288 single cells + 367 PICs
% 6. 'clustering13color' - clustering number for 6288 single cells
% 7. 'clustering13name' - cell typ names for 6288 single cells
% 8. 'clustering13name_unique' - 13 cell types
% 9. 'singletsIndex' - 6655 Boolean memberships for single cells
% 10. 'doubletsIndex' - 6655 Boolean memberships for PICs
% 11. 'sample_type_color' - sample type number for 6288 single cells + 367 PICs
%       (1 - E7.5 singlets, 2 - E8.5 singlets, 3 - E9.5 singlets,
%           4 - E7.5 doublets, 5 - E8.5 doublets, 6 - E9.5 doublets
clusterSize=max(clustering13color);
log_data_doublets=log_data_total(:,doubletsIndex);
%%%%%%%%%%%%% ---- Pipeline for clustering13 (E8.5) ---- %%%%%%%%%%%%%
folderName='mouseEmbryoE8.5';
log_data=log_data(:,sample_type_color(singletsIndex)==5);
clustering13color=clustering13color(sample_type_color(singletsIndex)==5);
log_data_doublets=log_data_doublets(:,sample_type_color(doubletsIndex)==2);
%%%%%%%%%%%%% Calculating z-value %%%%%%%%%%%%%
log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
log_data_doublets_zvalue=(log_data_doublets-repmat(mean(log_data_doublets,2),1,size(log_data_doublets,2)))./repmat(std(log_data_doublets')',1,size(log_data_doublets,2));

log_data_zvalue(isnan(log_data_zvalue))=0;
log_data_doublets_zvalue(isnan(log_data_doublets_zvalue))=0;
%%%%%%%%%%%%% 1) DEG by clustering 13 %%%%%%%%%%%%%

[pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data,clustering13color);
save([folderName,'/pvalue_fdr_logRatio_zvalue.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')

%%%%%%%%%%%%% 2) PIC SVM classification %%%%%%%%%%%%%
clusterSelect=[1,2,3,5,6,9,10];
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
pCutoff=0.01;lrCutoff=0.4;
seedNumber=1;randSize=5000;
DEGnumber=5;
[bestMatch,artificialDoubletsCombiUnique,SVMcl]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering13color,clusterSelect,clustering13name_unique,log_data_zvalue,log_data_doublets_zvalue,DEGnumber);
CV_SVMcl=crossval(SVMcl);
genError = kfoldLoss(CV_SVMcl)
save([folderName,'/SVM_bestMatch.mat'],'bestMatch','artificialDoubletsCombiUnique','SVMcl','CV_SVMcl','genError')

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

