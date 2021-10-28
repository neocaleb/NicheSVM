function [bestMatch,artificialDoubletsCombiUnique,SVMcl]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering,clusterSelect,clustering_name_unique,log_data_zvalue,log_data_doublets_zvalue,DEGnumber)

%%%%%%%%%%%%% 1) make artificial doublets %%%%%%%%%%%%%
clusterSize=size(pvalue_total,1);
DEGindex=zeros(size(log_data_zvalue,1),clusterSize);
for clusterIndex=1:clusterSize
    DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
end

geneIndex2=[];
for clusterIndex=1:size(clusterSelect,2)
    clusterIndex2=clusterSelect(clusterIndex);
    geneIndexTemp=find(DEGindex(:,clusterIndex2)&sum(DEGindex(:,clusterSelect),2)==1);
    [~,sortIndex]=sort(logRatio_total{clusterIndex2}(geneIndexTemp),'descend');
    if size(geneIndexTemp,1)>DEGnumber
        geneIndexTemp2=geneIndexTemp(sortIndex(1:5));
    else
        geneIndexTemp2=geneIndexTemp(sortIndex);
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
