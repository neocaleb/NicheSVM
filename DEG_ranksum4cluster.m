function [pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data,clustering)

pvalue_total=cell(clusterSize,1);
fdr_total=cell(clusterSize,1);
logRatio_total=cell(clusterSize,1);
zvalue_total=cell(clusterSize,1);
for clusterIndex=1:clusterSize
    pvalue=ones(size(log_data,1),1);
    fdr=ones(size(log_data,1),1);
    logRatio=zeros(size(log_data,1),1);
    zvalue=zeros(size(log_data,1),1);
    cellSelect1=clustering==clusterIndex;
    cellSelect2=clustering~=clusterIndex;
    if sum(cellSelect1)>0
        for i=1:size(log_data,1)
            [p,~,stats]=ranksum(log_data(i,find(cellSelect1)),log_data(i,find(cellSelect2)));
            pvalue(i)=p;
            logRatio(i)=mean(log_data(i,cellSelect1)+1)-mean(log_data(i,cellSelect2)+1);
            zvalue(i)=stats.zval;
        end
        fdr=mafdr(pvalue);
    end
    pvalue_total{clusterIndex}=pvalue;
    fdr_total{clusterIndex}=fdr;
    logRatio_total{clusterIndex}=logRatio;
    zvalue_total{clusterIndex}=zvalue;
end
