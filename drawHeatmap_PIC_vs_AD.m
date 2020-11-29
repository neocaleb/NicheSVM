function drawHeatmap_PIC_vs_AD(outputFolder,bestMatch,artificialDoubletsCombiUnique,clustering,clusterSelect,clustering_name_unique,log_data_zvalue,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name,DEGnumber,DEGnumberPIC_AD)

bestMatchUnique=unique(bestMatch);
bestMatchCount=histc(bestMatch,bestMatchUnique);
[~,bestMatchCountSortIndex]=sort(bestMatchCount,'descend');

clusterSize=size(pvalue_total,1);
DEGindex=zeros(size(gene_name,1),clusterSize);
for clusterIndex=1:clusterSize
    DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
end

bestMatchSize=size(bestMatchUnique,2);
for bestMatchIndex=1:bestMatchSize
    bestMatchUniqueTemp=bestMatchUnique(bestMatchCountSortIndex(bestMatchIndex));
    combiTemp=artificialDoubletsCombiUnique(bestMatchUniqueTemp).split('+');
    heteroIndex=find(bestMatch==bestMatchUniqueTemp);
    if size(heteroIndex,2)>=5 && ~strcmp(combiTemp{1},combiTemp{2})
        clusterIndex1=find(strcmp(clustering_name_unique,combiTemp(1)));
        clusterIndex2=find(strcmp(clustering_name_unique,combiTemp(2)));
        cellSelect1=find(clustering==clusterIndex1);
        cellSelect2=find(clustering==clusterIndex2);
        
        [~,geneIndex]=sort(logRatio_totalPIC_AD{bestMatchIndex},'descend');
        geneIndex=geneIndex(pvalue_totalPIC_AD{bestMatchIndex}(geneIndex)<pCutoff);
        if size(geneIndex,1)>DEGnumberPIC_AD
            geneIndex=geneIndex(1:DEGnumberPIC_AD);
        end
        
        [~,cellSortIndex]=sort(sum(log_data_doublets_zvalue(geneIndex,heteroIndex)),'descend');
        log_data_doublets_zvalueTemp=log_data_doublets_zvalue(:,heteroIndex);
        
        geneIndex1=find(DEGindex(:,clusterIndex1)&sum(DEGindex(:,clusterSelect),2)==1);
        [~,sortIndex1]=sort(logRatio_total{clusterIndex1}(geneIndex1),'descend');
        if size(sortIndex1,1)>DEGnumber
            geneIndex1=geneIndex1(sortIndex1(1:DEGnumber));
        else
            geneIndex1=geneIndex1(sortIndex1);
        end
        geneIndex2=find(DEGindex(:,clusterIndex2)&sum(DEGindex(:,clusterSelect),2)==1);
        [~,sortIndex2]=sort(logRatio_total{clusterIndex2}(geneIndex2),'descend');
        if size(sortIndex2,1)>DEGnumber
            geneIndex2=geneIndex2(sortIndex2(1:DEGnumber));
        else
            geneIndex2=geneIndex2(sortIndex2);
        end
        
        close all
        figure(1)
        subplot(3,3,1);
        imagesc(log_data_doublets_zvalueTemp(geneIndex,cellSortIndex))
        yticks(1:size(geneIndex,1))
        yticklabels(gene_name(geneIndex))
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        title(artificialDoubletsCombiUnique(bestMatchUniqueTemp),'FontSize',12)
        
        subplot(3,3,2);
        imagesc(log_data_zvalue(geneIndex,cellSelect1))
        yticks([])
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        title([clustering_name_unique{clusterIndex1},'(single cell)'],'FontSize',12)
        
        subplot(3,3,3);
        imagesc(log_data_zvalue(geneIndex,cellSelect2))
        yticks([])
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        title([clustering_name_unique{clusterIndex2},'(single cell)'],'FontSize',12)
        
        subplot(3,3,4);
        imagesc(log_data_doublets_zvalueTemp(geneIndex1,cellSortIndex))
        yticks(1:size(geneIndex1,1))
        yticklabels(gene_name(geneIndex1))
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        
        subplot(3,3,5);
        imagesc(log_data_zvalue(geneIndex1,cellSelect1))
        yticks([])
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        
        subplot(3,3,6);
        imagesc(log_data_zvalue(geneIndex1,cellSelect2))
        yticks([])
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        
        subplot(3,3,7);
        imagesc(log_data_doublets_zvalueTemp(geneIndex2,cellSortIndex))
        yticks(1:size(geneIndex2,1))
        yticklabels(gene_name(geneIndex2))
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        
        subplot(3,3,8);
        imagesc(log_data_zvalue(geneIndex2,cellSelect1))
        yticks([])
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        
        subplot(3,3,9);
        imagesc(log_data_zvalue(geneIndex2,cellSelect2))
        yticks([])
        xticks([])
        caxis([-3 3])
        colormap 'jet'
        
        set(gcf, 'Position', [100, 100, 550, 650])
        
        for i=7:9
            subplot(3,3,i)
            p = get(gca, 'Position');
            p(4) = p(4) / 1.2;
            set(gca, 'Position', p);
        end
        for i=4:6
            subplot(3,3,i)
            p = get(gca, 'Position');
            p(2) = p(2) - p(4)*0.5;
            p(4) = p(4) / 1.2;
            set(gca, 'Position', p);
        end
        for i=1:3
            subplot(3,3,i)
            p = get(gca, 'Position');
            p_diff = p(4) * 1;
            p(4) = p(4) + p_diff;
            p(2) = p(2) - p_diff;
            set(gca, 'Position', p);
        end
        
        
        saveas(gcf,[outputFolder,'/',artificialDoubletsCombiUnique{bestMatchUniqueTemp},'.pdf'])
    end
end
