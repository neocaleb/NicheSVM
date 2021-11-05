function DEGlists_PIC_vs_AD(outputFolder,bestMatch,artificialDoubletsCombiUnique,pCutoff,lrCutoff,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name)

bestMatchUnique=unique(bestMatch);
bestMatchCount=histc(bestMatch,bestMatchUnique);
[~,bestMatchCountSortIndex]=sort(bestMatchCount,'descend');

bestMatchSize=size(bestMatchUnique,2);
for bestMatchIndex=1:bestMatchSize
    bestMatchUniqueTemp=bestMatchUnique(bestMatchCountSortIndex(bestMatchIndex));
    combiTemp=artificialDoubletsCombiUnique(bestMatchUniqueTemp).split('+');
    heteroIndex=find(bestMatch==bestMatchUniqueTemp);
    if size(heteroIndex,2)>=5 && ~strcmp(combiTemp{1},combiTemp{2})
        logRatio=logRatio_totalPIC_AD{bestMatchIndex};
        pvalue=pvalue_totalPIC_AD{bestMatchIndex};
        
        [~,geneIndex]=sort(logRatio,'descend');
        geneIndex=geneIndex(pvalue(geneIndex)<pCutoff&logRatio(geneIndex)>lrCutoff);
        
        ofile = fopen([outputFolder,'/DEG_',artificialDoubletsCombiUnique{bestMatchUniqueTemp},'up_pvalue',num2str(pCutoff),'lr',num2str(lrCutoff),'.txt'],'w');
        ofile2 = fopen([outputFolder,'/DEG_',artificialDoubletsCombiUnique{bestMatchUniqueTemp},'up_pvalue',num2str(pCutoff),'lr',num2str(lrCutoff),'_with_pvalue_logratio.txt'],'w');
        for i=1:size(geneIndex,1)
            fprintf(ofile,'%s\n',gene_name(geneIndex(i)));
            fprintf(ofile2,'%s\t%f\t%f\n',gene_name(geneIndex(i)),pvalue(geneIndex(i)),abs(logRatio(geneIndex(i))));
        end
        fclose(ofile);
        fclose(ofile2);
    end
end
