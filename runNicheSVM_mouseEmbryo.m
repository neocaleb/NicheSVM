%%%%%%%%%%%%% Load MATLAB file %%%%%%%%%%%%%
load 'mouseEmbryoPIC.mat'
% 1. 'cell_id' - 6288 single cells
% 2. 'gene_name' - 1608 genes
% 3. 'log_data' - 1608 by 6288
% 4. 'clustering13color' - 6288 cells
% 5. 'clustering13name' - 6288 cells
% 6. 'clustering13name_unique' - 13 cell types
% 7. 'cell_id_total' - 6655 cells and PICs
% 8. 'log_data_total' - 6655 by 1608
% 9. 'singletsIndex' - 6655 Boolean memberships for single cells
% 10. 'doubletsIndex' - 6655 Boolean memberships for PICs
% 11. 'sample_type_color' - 6655 cells and PICs (1: PICs E7.5, 2: PICs for
% E8.5, 3: PICs for E9.5, 4: single cells for E7.5, 5: single cells for
% E8.5, 6: single cells for E9.5)

%%%%%%%%%%%%% Calculating z-value %%%%%%%%%%%%%
log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
log_data_doublets=log_data_total(:,doubletsIndex);
log_data_doublets_zvalue=(log_data_doublets-repmat(mean(log_data_doublets,2),1,size(log_data_doublets,2)))./repmat(std(log_data_doublets')',1,size(log_data_doublets,2));

%%%%%%%%%%%%% ---- Pipeline for clustering13 ---- %%%%%%%%%%%%%
folderName='mouseEmbryo';
%%%%%%%%%%%%% DEG by clustering 13 (totla and each time point) %%%%%%%%%%%%%
clusterSize=max(clustering13color);
cellSelect=sample_type_color(singletsIndex)>3; % total
[pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data(:,cellSelect),clustering13color(cellSelect));
save([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')
cellSelect=sample_type_color(singletsIndex)==4; % E7.5
[pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data(:,cellSelect),clustering13color(cellSelect));
save([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13E7_5.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')
cellSelect=sample_type_color(singletsIndex)==5; % E8.5
[pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data(:,cellSelect),clustering13color(cellSelect));
save([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13E8_5.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')
cellSelect=sample_type_color(singletsIndex)==6; % E9.5
[pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data(:,cellSelect),clustering13color(cellSelect));
save([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13E7_5.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')

%%%%%%%%%%%%% --- PIC SVM classification (E7.5) --- %%%%%%%%%%%%%
clusterSelect=[3,7,9,11,12,13];
cellSelect1=sample_type_color(singletsIndex)==4;
cellSelect2=sample_type_color(doubletsIndex)==1;
load([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13E7_5.mat'])
pCutoff=0.01;lrCutoff=0.4;
seedNumber=1;randSize=10000;
DEGnumber=5;
[bestMatch,artificialDoubletsCombiUnique]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering13color(cellSelect1),clusterSelect,clustering13name_unique,log_data_zvalue(:,cellSelect1),log_data_doublets_zvalue(:,cellSelect2),DEGnumber);
save([folderName,'/SVM_bestMatchE7_5.mat'],'bestMatch','artificialDoubletsCombiUnique')

%%%%%%%%%%%%% --- PIC SVM classification (E8.5) --- %%%%%%%%%%%%%
clusterSelect=[1,2,3,5,6,9,10];
cellSelect1=sample_type_color(singletsIndex)==5;
cellSelect2=sample_type_color(doubletsIndex)==2;
load([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13E8_5.mat'])
pCutoff=0.01;lrCutoff=0.4;
seedNumber=1;randSize=10000;
DEGnumber=5;
[bestMatch,artificialDoubletsCombiUnique]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering13color(cellSelect1),clusterSelect,clustering13name_unique,log_data_zvalue(:,cellSelect1),log_data_doublets_zvalue(:,cellSelect2),DEGnumber);
save([folderName,'/SVM_bestMatchE8_5.mat'],'bestMatch','artificialDoubletsCombiUnique')

%%%%%%%%%%%%% --- PIC SVM classification (E9.5) --- %%%%%%%%%%%%%
clusterSelect=[2,5,6,8,10];
cellSelect1=sample_type_color(singletsIndex)==6;
cellSelect2=sample_type_color(doubletsIndex)==3;
load([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13E9_5.mat'])
pCutoff=0.01;lrCutoff=0.4;
seedNumber=1;randSize=10000;
DEGnumber=5;
[bestMatch,artificialDoubletsCombiUnique]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering13color(cellSelect1),clusterSelect,clustering13name_unique,log_data_zvalue(:,cellSelect1),log_data_doublets_zvalue(:,cellSelect2),DEGnumber);
save([folderName,'/SVM_bestMatchE9_5.mat'],'bestMatch','artificialDoubletsCombiUnique')

%%%%%%%%%%%%% --- Integrating PIC SVM classification results (E7.5 + E8.5 + E9.5) --- %%%%%%%%%%%%%
clusterSelect=[1,2,3,5,6,7,8,9,10,11,12,13];
seedNumber=1;randSize=10000;
[~,~,artificialDoubletsCombiUniqueCombined,~]=generateAD(seedNumber,randSize,clustering13color,clusterSelect,log_data_zvalue(1,:),clustering13name_unique);

bestMatchCombined=[];
load([folderName,'/SVM_bestMatchE7_5.mat'])
for i=1:size(bestMatch,2)
    bestMatch(i)=find(strcmp(artificialDoubletsCombiUnique(bestMatch(i)),artificialDoubletsCombiUniqueCombined));
end
bestMatchCombined=[bestMatchCombined bestMatch];
load([folderName,'/SVM_bestMatchE8_5.mat'])
for i=1:size(bestMatch,2)
    bestMatch(i)=find(strcmp(artificialDoubletsCombiUnique(bestMatch(i)),artificialDoubletsCombiUniqueCombined));
end
bestMatchCombined=[bestMatchCombined bestMatch];
load([folderName,'/SVM_bestMatchE9_5.mat'])
for i=1:size(bestMatch,2)
    bestMatch(i)=find(strcmp(artificialDoubletsCombiUnique(bestMatch(i)),artificialDoubletsCombiUniqueCombined));
end
bestMatchCombined=[bestMatchCombined bestMatch];

artificialDoubletsCombiUnique=artificialDoubletsCombiUniqueCombined;
bestMatch=bestMatchCombined;
clear artificialDoubletsCombiUniqueCombined bestMatchCombined

%%% Draw heatmap
load([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13.mat'])
clusterSelect=[1,2,3,5,6,7,8,9,10,11,12,13];
outputFile=[folderName,'/heatmap_PICSVM_top5DEG.pdf'];
DEGnumber=5;
drawHeatmap_PICSVM(outputFile,bestMatch,artificialDoubletsCombiUnique,clusterSelect,clustering13name_unique,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,gene_name,DEGnumber);

%%%%%%%%%%%%% DEG between PICs and artificial doublets %%%%%%%%%%%%%
clusterSelect=[1,2,3,5,6,7,8,9,10,11,12,13];
seedNumber=1;randSize=10000;
minCell=5;
[pvalue_total,fdr_total,logRatio_total]=DEG_PIC_vs_AD(bestMatch,seedNumber,randSize,clustering13color,clusterSelect,log_data_zvalue,clustering13name_unique,log_data_doublets_zvalue,minCell);
save([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'],'pvalue_total','fdr_total','logRatio_total')

outputFolder=[folderName,'/heatmapPIC_vs_AD'];
load([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'])
pvalue_totalPIC_AD=pvalue_total;
logRatio_totalPIC_AD=logRatio_total;
load([folderName,'/pvalue_fdr_logRatio_zvalue_clustering13.mat'])
DEGnumber=5;DEGnumberPIC_AD=10;
drawHeatmap_PIC_vs_AD(outputFolder,bestMatch,artificialDoubletsCombiUnique,clustering13color,clusterSelect,clustering13name_unique,log_data_zvalue,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name,DEGnumber,DEGnumberPIC_AD);

