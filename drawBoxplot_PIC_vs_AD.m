function drawBoxplot_PIC_vs_AD(seedNumber,randSize,outputFolder,bestMatch,artificialDoubletsCombiUnique,clustering,clusterSelect,clustering_name_unique,log_data_zvalue,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name,DEGnumber,DEGnumberPIC_AD)

bestMatchUnique=unique(bestMatch);
bestMatchCount=histc(bestMatch,bestMatchUnique);
[~,bestMatchCountSortIndex]=sort(bestMatchCount,'descend');

clusterSize=size(pvalue_total,1);
DEGindex=zeros(size(gene_name,1),clusterSize);
for clusterIndex=1:clusterSize
    DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
end

[~,artificialDoubletsCombiColor,~,~]=generateAD(seedNumber,randSize,clustering,clusterSelect,log_data_zvalue(1,:),clustering_name_unique);

bestMatchSize=size(bestMatchUnique,2);
for bestMatchIndex=1:bestMatchSize
    bestMatchUniqueTemp=bestMatchUnique(bestMatchCountSortIndex(bestMatchIndex));
    combiTemp=artificialDoubletsCombiUnique(bestMatchUniqueTemp).split('+');
    heteroIndex=find(bestMatch==bestMatchUniqueTemp);
    if size(heteroIndex,2)>=5 && ~strcmp(combiTemp{1},combiTemp{2})
        clusterIndex1=find(strcmp(clustering_name_unique,combiTemp(1)));
        clusterIndex2=find(strcmp(clustering_name_unique,combiTemp(2)));
        
        [~,geneIndex]=sort(logRatio_totalPIC_AD{bestMatchIndex},'descend');
        geneIndex=geneIndex(pvalue_totalPIC_AD{bestMatchIndex}(geneIndex)<pCutoff);
        if size(geneIndex,1)>DEGnumberPIC_AD
            geneIndex=geneIndex(1:DEGnumberPIC_AD);
        end
        
        geneIndexTotal=cell(2,1);
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
        geneIndexTotal{1}=geneIndex1;
        geneIndexTotal{2}=geneIndex2;
        
        close all
        
        figure(1)
        for j=1:DEGnumberPIC_AD
            [log_data_artificialDoublets,~,~,~]=generateAD(seedNumber,randSize,clustering,clusterSelect,log_data_zvalue(geneIndex(j),:),clustering_name_unique);
        
            boxplotData=[log_data_doublets_zvalue(geneIndex(j),heteroIndex),log_data_artificialDoublets(:,artificialDoubletsCombiColor==bestMatchUniqueTemp)];
            boxplotDataLabel=strings(size(boxplotData));
            boxplotDataLabel(1:size(heteroIndex,2))='Observed';
            boxplotDataLabel(size(heteroIndex,2)+1:end)='Expected';
            boxSize1=quantile(boxplotData(1:size(heteroIndex,2)),0.75)-quantile(boxplotData(1:size(heteroIndex,2)),0.25);
            whiskerUpper1=quantile(boxplotData(1:size(heteroIndex,2)),0.75)+boxSize1*2;
            whiskerLower1=quantile(boxplotData(1:size(heteroIndex,2)),0.25)-boxSize1*2;
            boxSize2=quantile(boxplotData(size(heteroIndex,2)+1:end),0.75)-quantile(boxplotData(size(heteroIndex,2)+1:end),0.25);
            whiskerUpper2=quantile(boxplotData(size(heteroIndex,2)+1:end),0.75)+boxSize2*2;
            whiskerLower2=quantile(boxplotData(size(heteroIndex,2)+1:end),0.25)-boxSize2*2;
            subplot(2,DEGnumberPIC_AD/2,j)
            boxplot(boxplotData,boxplotDataLabel,'Symbol','')
            ylim([min(whiskerLower1,whiskerLower2) max(whiskerUpper1,whiskerUpper2)])
            xticks([])
            title(gene_name(geneIndex(j)),'FontSize',12)
        end
        set(gcf, 'Position', [100, 100, 600, 180])
        saveas(gcf,[outputFolder,'/',artificialDoubletsCombiUnique{bestMatchUniqueTemp},'.NeighborMarker.boxplot.pdf'])
        
        figure(2)
        for geneIndexTempIndex=1:2
            geneIndexTemp=geneIndexTotal{geneIndexTempIndex};
            for j=1:DEGnumber
                if j<=size(geneIndexTemp,1)
                    [log_data_artificialDoublets,~,~,~]=generateAD(seedNumber,randSize,clustering,clusterSelect,log_data_zvalue(geneIndexTemp(j),:),clustering_name_unique);
                    boxplotData=[log_data_doublets_zvalue(geneIndexTemp(j),heteroIndex),log_data_artificialDoublets(:,artificialDoubletsCombiColor==bestMatchUniqueTemp)];
                    boxplotDataLabel=strings(size(boxplotData));
                    boxplotDataLabel(1:size(heteroIndex,2))='Observed';
                    boxplotDataLabel(size(heteroIndex,2)+1:end)='Expected';
                    boxSize1=quantile(boxplotData(1:size(heteroIndex,2)),0.75)-quantile(boxplotData(1:size(heteroIndex,2)),0.25);
                    whiskerUpper1=quantile(boxplotData(1:size(heteroIndex,2)),0.75)+boxSize1*2;
                    whiskerLower1=quantile(boxplotData(1:size(heteroIndex,2)),0.25)-boxSize1*2;
                    boxSize2=quantile(boxplotData(size(heteroIndex,2)+1:end),0.75)-quantile(boxplotData(size(heteroIndex,2)+1:end),0.25);
                    whiskerUpper2=quantile(boxplotData(size(heteroIndex,2)+1:end),0.75)+boxSize2*2;
                    whiskerLower2=quantile(boxplotData(size(heteroIndex,2)+1:end),0.25)-boxSize2*2;
                    subplot(2,DEGnumber,j+(geneIndexTempIndex-1)*DEGnumber)
                    boxplot(boxplotData,boxplotDataLabel,'symbol','')
                    ylim([min(whiskerLower1,whiskerLower2) max(whiskerUpper1,whiskerUpper2)])
                    title(gene_name(geneIndexTemp(j)),'FontSize',12)
                    xticks([])
                end
            end
        end
        set(gcf, 'Position', [100, 100, 600, 180])
        saveas(gcf,[outputFolder,'/',artificialDoubletsCombiUnique{bestMatchUniqueTemp},'.CellTypeMarker.boxplot.pdf'])
        
    end
end
