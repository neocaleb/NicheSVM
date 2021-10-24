%%%%%%%%%%%%% Load MATLAB file %%%%%%%%%%%%%
load 'liver_paired.mat'
% 1. 'cell_id' - 3719 single cells
% 2. 'cell_id_total' - 3719 single cells + 4621 PICs
% 3. 'gene_name_total' - 33948 genes
% 4. 'log_data_total' - 33948 by 3719 single cells + 4621 PICs
% 5. 'clustering7color' - clustering number for 3719 single cells
% 6. 'clustering7name' - cell type names for 3719 single cells
% 7. 'clustering7name_unique' - 7 cell types
% 8. 'singletsIndex' - 8340 Boolean memberships for single cells
% 9. 'doubletsIndex' - 8340 Boolean memberships for PICs
gene_name=gene_name_total;
%%%%%%%%%%%%% Calculating z-value %%%%%%%%%%%%%
log_data=log_data_total(:,singletsIndex);

log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
log_data_doublets=log_data_total(:,doubletsIndex);
log_data_doublets_zvalue=(log_data_doublets-repmat(mean(log_data_doublets,2),1,size(log_data_doublets,2)))./repmat(std(log_data_doublets')',1,size(log_data_doublets,2));

log_data_zvalue(isnan(log_data_zvalue))=0;
log_data_doublets_zvalue(isnan(log_data_doublets_zvalue))=0;

%%%%%%%%%%%%% ---- Pipeline for clustering7 ---- %%%%%%%%%%%%%
folderName='liverPaired';
%%%%%%%%%%%%% 1) DEG by clustering 7 %%%%%%%%%%%%%
clusterSize=max(clustering7color);
[pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data,clustering7color);
save([folderName,'/pvalue_fdr_logRatio_zvalue.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')

%%%%%%%%%%%%% 2) PIC SVM classification %%%%%%%%%%%%%
clusterSelect=1:7;
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
pCutoff=0.01;lrCutoff=0.5;
DEGindex=zeros(size(gene_name,1),clusterSize);
for clusterIndex=1:clusterSize
    cellIndexTemp=find(clustering7color==clusterIndex);
    DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
end
sum(DEGindex)

seedNumber=1;randSize=10000;
DEGnumber=5;
[bestMatch,artificialDoubletsCombiUnique]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering7color,clusterSelect,clustering7name_unique,log_data_zvalue,log_data_doublets_zvalue,DEGnumber);
save([folderName,'/SVM_bestMatch.mat'],'bestMatch','artificialDoubletsCombiUnique')

%%% Draw heatmap
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
load([folderName,'/SVM_bestMatch.mat'])
clusterSelect=1:7;
outputFile=[folderName,'/heatmap_PICSVM_top5DEG.pdf'];
pCutoff=0.01;lrCutoff=0.5;
DEGnumber=5;
drawHeatmap_PICSVM(outputFile,bestMatch,artificialDoubletsCombiUnique,clusterSelect,clustering7name_unique,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,gene_name,DEGnumber);

%%%%%%%%%%%%% 3) DEG between PICs and artificial doublets %%%%%%%%%%%%%
clusterSelect=1:7;
seedNumber=1;randSize=10000;
minCell=5;
[pvalue_total,fdr_total,logRatio_total]=DEG_PIC_vs_AD(bestMatch,seedNumber,randSize,clustering7color,clusterSelect,log_data_zvalue,clustering7name_unique,log_data_doublets_zvalue,minCell);
save([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'],'pvalue_total','fdr_total','logRatio_total')

outputFolder=folderName;
load([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'])
pvalue_totalPIC_AD=pvalue_total;
logRatio_totalPIC_AD=logRatio_total;
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
pCutoff=0.01;lrCutoff=0.5;
DEGnumber=5;DEGnumberPIC_AD=10;
drawHeatmap_PIC_vs_AD(outputFolder,bestMatch,artificialDoubletsCombiUnique,clustering7color,clusterSelect,clustering7name_unique,log_data_zvalue,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name,DEGnumber,DEGnumberPIC_AD);

