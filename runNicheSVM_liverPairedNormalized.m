%%%%%%%%%%%%% Load MATLAB file %%%%%%%%%%%%%
load 'liver_paired_normalized.mat'
% 1. 'cell_id' - 4568 single cells
% 2. 'cell_id_total' - 4568 single cells + 4621 PICs
% 3. 'gene_name_filtered' - 2999 genes
% 4. 'log_data_total_normalized' - 2999 by 4568 single cells + 4621 PICs
% 5. 'clustering8color' - clustering number for 4568 single cells
% 6. 'clustering8name' - cell type names for 4568 single cells
% 7. 'clustering8name_unique' - 8 cell types
% 8. 'singletsIndex' - 9189 Boolean memberships for single cells
% 9. 'doubletsIndex' - 9189 Boolean memberships for PICs
gene_name=gene_name_filtered;
log_data_total=log_data_total_normalized;

%%%%%%%%%%%%% gene filter %%%%%%%%%%%%%
ERCCindex=zeros(size(gene_name));
for i=1:size(gene_name,1)
    if size(gene_name{i},2)>3
        if sum(gene_name{i}(1:4)=='ERCC')==4
            ERCCindex(i)=1;
        end
    end
end
ERCCindex=ERCCindex==1;

existCutoff=0.01;
existIndex=sum(log_data_total(:,singletsIndex)>1,2)>sum(singletsIndex)*existCutoff&sum(log_data_total(:,doubletsIndex)>1,2)>sum(doubletsIndex)*existCutoff;

log_data_total=log_data_total(existIndex&~ERCCindex,:);
gene_name=gene_name(existIndex&~ERCCindex);

%%%%%%%%%%%%% Calculating z-value %%%%%%%%%%%%%
log_data=log_data_total(:,singletsIndex);

log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
log_data_doublets=log_data_total(:,doubletsIndex);
log_data_doublets_zvalue=(log_data_doublets-repmat(mean(log_data_doublets,2),1,size(log_data_doublets,2)))./repmat(std(log_data_doublets')',1,size(log_data_doublets,2));

log_data_zvalue(isnan(log_data_zvalue))=0;
log_data_doublets_zvalue(isnan(log_data_doublets_zvalue))=0;

%%%%%%%%%%%%% ---- Pipeline for clustering8 ---- %%%%%%%%%%%%%
folderName='liverPairedNormalized';
%%%%%%%%%%%%% 1) DEG by clustering 8 %%%%%%%%%%%%%
clusterSize=max(clustering8color);
[pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data,clustering8color);
save([folderName,'/pvalue_fdr_logRatio_zvalue.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')

%%%%%%%%%%%%% 2) PIC SVM classification %%%%%%%%%%%%%
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
pCutoff=0.01;lrCutoff=0.3;
DEGindex=zeros(size(gene_name,1),clusterSize);
for clusterIndex=1:clusterSize
    cellIndexTemp=find(clustering6color==clusterIndex);
    DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
end
DEGindexOnly=zeros(size(gene_name,1),clusterSize);
clusterSelect=1:clusterSize;
for clusterIndex=1:size(clusterSelect,2)
    clusterIndex=clusterSelect(clusterIndex);
    for i=1:size(gene_name,1)
        DEGindexOnly(i,clusterIndex)=DEGindex(i,clusterIndex) && sum(DEGindex(i,clusterSelect),2)==1;
        if DEGindex(i,clusterIndex) && sum(DEGindex(i,clusterSelect),2)==2
            clusterTemp=clusterSelect(find(DEGindex(i,clusterSelect)));
            DEGindexOnly(i,clusterIndex)=logRatio_total{clusterIndex}(i)-logRatio_total{clusterTemp(clusterTemp~=clusterIndex)}(i)>lrCutoff;
        end
    end
end
sum(DEGindexOnly)

clusterSelect=find(sum(DEGindexOnly)>1);

seedNumber=1;randSize=10000;
DEGnumber=5;
[bestMatch,artificialDoubletsCombiUnique,SVMcl]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering8color,clusterSelect,clustering8name_unique,log_data_zvalue,log_data_doublets_zvalue,DEGnumber);
save([folderName,'/SVM_bestMatch.mat'],'SVMcl','bestMatch','artificialDoubletsCombiUnique')

%%% Draw heatmap
load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
load([folderName,'/SVM_bestMatch.mat'])
outputFile=[folderName,'/heatmap_PICSVM_top5DEG.pdf'];
pCutoff=0.01;lrCutoff=0.3;
DEGnumber=10;
drawHeatmap_PICSVM(outputFile,bestMatch,artificialDoubletsCombiUnique,clusterSelect,clustering8name_unique,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,gene_name,DEGnumber);

%%%%%%%%%%%%% 3) DEG between PICs and artificial doublets %%%%%%%%%%%%%
load([folderName,'/SVM_bestMatch.mat'])
clusterSelect=find(sum(DEGindexOnly)>1);
seedNumber=1;randSize=10000;
minCell=5;
[pvalue_total,fdr_total,logRatio_total]=DEG_PIC_vs_AD(bestMatch,seedNumber,randSize,clustering8color,clusterSelect,log_data_zvalue,clustering8name_unique,log_data_doublets_zvalue,minCell);
save([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'],'pvalue_total','fdr_total','logRatio_total')

outputFolder=folderName;
load([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'])
pvalue_totalPIC_AD=pvalue_total;
logRatio_totalPIC_AD=logRatio_total;
load([folderName2,'/pvalue_fdr_logRatio_zvalue.mat'])
pCutoff=0.01;lrCutoff=0.5;
%%% Save DEG lists
DEGlists_PIC_vs_AD(outputFolder,bestMatch,artificialDoubletsCombiUnique,pCutoff,lrCutoff,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name);
%%% boxplots for observed and expected expression of the neighbor-specific markers and cell type markers
DEGnumber=5;DEGnumberPIC_AD=10;
drawBoxplot_PIC_vs_AD(seedNumber,randSize,outputFolder,bestMatch,artificialDoubletsCombiUnique,clustering8color,clusterSelect,clustering8name_unique,log_data_zvalue,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name,DEGnumber,DEGnumberPIC_AD);
%%% heatmaps for the neighbor-specific markers and cell type markers
DEGnumber=5;DEGnumberPIC_AD=10;
drawHeatmap_PIC_vs_AD(outputFolder,bestMatch,artificialDoubletsCombiUnique,clustering8color,clusterSelect,clustering8name_unique,log_data_zvalue,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name,DEGnumber,DEGnumberPIC_AD);

