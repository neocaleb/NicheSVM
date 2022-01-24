function [bestMatch,artificialDoubletsCombiUnique,SVMcl]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering,clusterSelect,clustering_name_unique,log_data_zvalue,log_data_doublets_zvalue,DEGnumber)

%%%%%%%%%%%%% 1) make artificial doublets %%%%%%%%%%%%%
clusterSize=size(pvalue_total,1);
DEGindex=zeros(size(log_data_zvalue,1),clusterSize);
for clusterIndex=1:clusterSize
    DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
end

DEGindexOnly=zeros(size(gene_name,1),clusterSize);
for clusterIndex=1:size(clusterSelect,2)
    for i=1:size(gene_name,1)
        DEGindexOnly(i,clusterIndex)=DEGindex(i,clusterIndex) && sum(DEGindex(i,clusterSelect),2)==1;
        if DEGindex(i,clusterIndex) && sum(DEGindex(i,clusterSelect),2)==2
            clusterTemp=clusterSelect(find(DEGindex(i,clusterSelect)));
            DEGindexOnly(i,clusterIndex)=logRatio_total{clusterIndex}(i)-logRatio_total{clusterTemp(clusterTemp~=clusterIndex)}(i)>lrCutoff;
        end
    end
end

geneIndex=[];geneIndex2=[];
for clusterIndex=1:size(clusterSelect,2)
    clusterIndex=clusterSelect(clusterIndex);
    geneIndexTemp=find(DEGindexOnly(:,clusterIndex));
    [~,sortIndex]=sort(logRatio_total{clusterIndex}(geneIndexTemp),'descend');
    geneIndexTemp=geneIndexTemp(sortIndex);
    geneIndex=[geneIndex;geneIndexTemp];
    if size(geneIndexTemp,1)>DEGnumber
        geneIndexTemp2=geneIndexTemp(1:DEGnumber);
    else
        geneIndexTemp2=geneIndexTemp;
    end
    geneIndex2=[geneIndex2;geneIndexTemp2];
end


[log_data_artificialDoublets,artificialDoubletsCombiColor,artificialDoubletsCombiUnique,~]=generateAD(seedNumber,randSize,clustering,clusterSelect,log_data_zvalue(geneIndex2,:),clustering_name_unique);

%%%%%%%%%%%%% 2) run SVM %%%%%%%%%%%%%
SVMcl = fitcecoc(log_data_artificialDoublets',artificialDoubletsCombiColor');
[SVMlabels,~] = predict(SVMcl,log_data_doublets_zvalue(geneIndex2,:)');

bestMatch=SVMlabels';
bestMatchUnique=unique(bestMatch);
bestMatchCount=histc(bestMatch,bestMatchUnique);
[~,bestMatchCountSortIndex]=sort(bestMatchCount,'descend');

[["PIC type";artificialDoubletsCombiUnique(bestMatchUnique(bestMatchCountSortIndex))'],["Number";bestMatchCount(bestMatchCountSortIndex)']]
