function [log_data_artificialDoublets,artificialDoubletsCombiColor,artificialDoubletsCombiUnique,artificialDoubletsCombiUniqueLabelIndex]=generateAD(seedNumber,randSize,clustering,clusterSelect,log_data,clustering_name_unique)

rand('seed',seedNumber);
clusterSize=size(clustering_name_unique,2);clusterSize2=size(clusterSelect,2);
artificialSize=clusterSize2*(clusterSize2-1)/2+clusterSize2;
log_data_artificialDoublets=zeros(size(log_data,1),artificialSize*randSize);
artificialDoubletsCombiColor=zeros(1,artificialSize*randSize);
artificialDoubletsCombiUnique=strings(1,artificialSize);
artificialDoubletsCombiUniqueLabelIndex=zeros(1,artificialSize);
doubletsIteration=1;
for clusterIndex1=1:clusterSize
    for clusterIndex2=1:clusterSize
        if clusterIndex1<=clusterIndex2 && sum(clusterSelect==clusterIndex1)>0 && sum(clusterSelect==clusterIndex2)>0
            cellIndex1=find(clustering==clusterIndex1);
            cellIndex2=find(clustering==clusterIndex2);
            randIndex1=randi([1 size(cellIndex1,2)],1,randSize);
            randIndex2=randi([1 size(cellIndex2,2)],1,randSize);
            data_temp=(log_data(:,cellIndex1(randIndex1))+log_data(:,cellIndex2(randIndex2)))/2;
            log_data_artificialDoublets(:,(doubletsIteration-1)*randSize+1:doubletsIteration*randSize)=data_temp;
            artificialDoubletsCombiColor((doubletsIteration-1)*randSize+1:doubletsIteration*randSize)=doubletsIteration;
            artificialDoubletsCombiUnique(doubletsIteration)=[clustering_name_unique{clusterIndex1},'+',clustering_name_unique{clusterIndex2}];
            artificialDoubletsCombiUniqueLabelIndex(doubletsIteration)=randSize*(doubletsIteration-1)+round(randSize/2);
            doubletsIteration=doubletsIteration+1;
        end
    end
end