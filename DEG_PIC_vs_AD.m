function [pvalue_total,fdr_total,logRatio_total]=DEG_PIC_vs_AD(bestMatch,seedNumber,randSize,clustering,clusterSelect,log_data_zvalue,clustering_name_unique,log_data_doublets_zvalue,minCell)

bestMatchUnique=unique(bestMatch);
bestMatchCount=histc(bestMatch,bestMatchUnique);
[~,bestMatchCountSortIndex]=sort(bestMatchCount,'descend');

[~,artificialDoubletsCombiColor,~,~]=generateAD(seedNumber,randSize,clustering,clusterSelect,log_data_zvalue(1,:),clustering_name_unique);

bestMatchSize=size(bestMatchUnique,2);
pvalue_total=cell(bestMatchSize,1);
fdr_total=cell(bestMatchSize,1);
logRatio_total=cell(bestMatchSize,1);
for bestMatchIndex=1:bestMatchSize
    bestMatchUniqueTemp=bestMatchUnique(bestMatchCountSortIndex(bestMatchIndex));
    heteroIndex=find(bestMatch==bestMatchUniqueTemp);
    pvalue=ones(size(log_data_zvalue,1),1);
    fdr=ones(size(log_data_zvalue,1),1);
    logRatio=zeros(size(log_data_zvalue,1),1);
    if size(heteroIndex,2)>=minCell
        for i=1:size(log_data_zvalue,1)
            [log_data_artificialDoublets,~,~,~]=generateAD(seedNumber,randSize,clustering,clusterSelect,log_data_zvalue(i,:),clustering_name_unique);
            [p,h,stats]=ranksum(log_data_doublets_zvalue(i,heteroIndex),log_data_artificialDoublets(:,artificialDoubletsCombiColor==bestMatchUniqueTemp));
            pvalue(i,1)=p;
            logRatio(i,1)=mean(log_data_doublets_zvalue(i,heteroIndex))-mean(log_data_artificialDoublets(:,artificialDoubletsCombiColor==bestMatchUniqueTemp));
        end
        fdr=mafdr(pvalue);
    end
        
    pvalue_total{bestMatchIndex}=pvalue;
    fdr_total{bestMatchIndex}=fdr;
    logRatio_total{bestMatchIndex}=logRatio;
end