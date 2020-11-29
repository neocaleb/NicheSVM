function drawHeatmap_PICSVM(outputFile,bestMatch,artificialDoubletsCombiUnique,clusterSelect,clustering_name_unique,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,gene_name,DEGnumber)

load 'colormap_2to19grey.mat'

bestMatchUnique=unique(bestMatch);
bestMatchCount=histc(bestMatch,bestMatchUnique);
[~,bestMatchCountSortIndex]=sort(bestMatchCount,'descend');
bestMatchSortIndex=[];
bestMatchLabels=[];
bestMatchLabelIndex=[];
bestMatchSingleColor=zeros(2,size(bestMatch,2));
for i=1:size(bestMatchCountSortIndex,2)
    indexTemp=find(bestMatch==bestMatchUnique(bestMatchCountSortIndex(i)));
    bestMatchLabels=[bestMatchLabels artificialDoubletsCombiUnique(bestMatchUnique(bestMatchCountSortIndex(i)))];
    bestMatchLabelIndex=[bestMatchLabelIndex size(bestMatchSortIndex,2)+round(size(indexTemp,2)/2)];
        
    bestMatchSortIndex=[bestMatchSortIndex indexTemp];
    combiTemp=artificialDoubletsCombiUnique(bestMatchUnique(bestMatchCountSortIndex(i))).split('+');
    clusterIndex1=find(strcmp(clustering_name_unique,combiTemp(1)));
    clusterIndex2=find(strcmp(clustering_name_unique,combiTemp(2)));
    bestMatchSingleColor(1,indexTemp)=clusterIndex1;
    bestMatchSingleColor(2,indexTemp)=clusterIndex2;
end

clusterSize=size(pvalue_total,1);
DEGindex=zeros(size(gene_name,1),clusterSize);
for clusterIndex=1:clusterSize
    DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
end

geneIndex2=[];
for clusterIndex=1:size(clusterSelect,2)
    clusterIndex2=clusterSelect(clusterIndex);
    geneIndexTemp=find(DEGindex(:,clusterIndex2)&sum(DEGindex(:,clusterSelect),2)==1);
    [~,sortIndex]=sort(logRatio_total{clusterIndex2}(geneIndexTemp),'descend');
    if size(geneIndexTemp,1)>DEGnumber
        geneIndexTemp2=geneIndexTemp(sortIndex(1:DEGnumber));
    else
        geneIndexTemp2=geneIndexTemp(sortIndex);
    end
    geneIndex2=[geneIndex2;geneIndexTemp2];
end

log_data_temp=log_data_doublets_zvalue(geneIndex2,:);
close all
figure(1)
ax(1)=subplot(2,1,1);
imagesc(log_data_temp(:,bestMatchSortIndex))
caxis([-3 3])
colormap(ax(1),'jet')
yticks([1:size(geneIndex2,1)])
yticklabels(gene_name(geneIndex2));
xticks([])
set(gcf, 'Position', [100, 100, 600, 500])
ax(2)=subplot(2,1,2);
imagesc(bestMatchSingleColor(:,bestMatchSortIndex))
caxis([0 13])
colormap(ax(2),mycmap2to19grey{12})
yticks([])
xticks(bestMatchLabelIndex)
xticklabels(bestMatchLabels)
xtickangle(90)
subplot(2,1,2)
p = get(gca, 'Position');
p(4) = p(4) / 8;
set(gca, 'Position', p);
subplot(2,1,1)
p = get(gca, 'Position');
p_diff = p(4) * 1.2;
p(4) = p(4) + p_diff;
p(2) = p(2) - p_diff;
set(gca, 'Position', p);
saveas(gcf,outputFile)

