function TrainSet = TrainSetUpdate(TrainSet,Archive,Population,Offspring)
% Updating the training set of EMMA-GNG.

%------------------------------- Reference --------------------------------
% Liu Y, Zhang L, Zeng X, et al. Evolutionary multimodal multiobjective 
% optimization guided by growing neural gas[J]. Swarm and Evolutionary 
% Computation, 2024, 86: 101500.
%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    TrainSet = [TrainSet,Population,Offspring];
    Com = [Archive,TrainSet];
    [frontNo,maxNo] = NDSort(Com.objs,1);
    indexSet = setdiff(find(frontNo==1),1:size(Archive,2));
    Set = Com(indexSet);     

    Dim = size(TrainSet.objs,2);
    if Dim == 2
        thresholdObj = 1*10^(-3);
    end
    if Dim ==3
        thresholdObj = 1*10^(-2);
    end
    if Dim > 3
        thresholdObj = 0.2;
    end
    temp = [Population,Offspring];
    DataDisObj = pdist2(Archive.objs,temp.objs);
    [MinDataDisObj,~] = min(DataDisObj,[],1);
    Index = MinDataDisObj <= thresholdObj;
    temp = temp(Index);
    TrainSet = [temp,Set];
    [~,ia] = unique(TrainSet.decs,'rows','stable');
    TrainSet = TrainSet(ia);
 
end
