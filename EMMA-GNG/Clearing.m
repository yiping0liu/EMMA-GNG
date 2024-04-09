function TrainSet = Clearing(TrainSet) 
    [No,maxNo] = NDSort(TrainSet.objs,inf);
    [~,index] = sort(No,'ascend');
    TrainSet = TrainSet(index);
    Choose = true(1,length(TrainSet));
    sigma_dec = 0.002;
    dim = size(TrainSet.objs,2);
    
    if dim ==3
        sigma_dec = 0.02;
    elseif dim > 3
        sigma_dec = 0.04;
    end

    Select = Clear(TrainSet,sigma_dec,Choose);   
    Choose = Select;
    TrainSet = TrainSet(Choose);

end
function Choose = Clear(TrainSet,sigma_radius,Choose)
    RemoveSet = [];
    dis = pdist2(TrainSet.decs,TrainSet.decs);
    for i = 1:length(TrainSet)-1
        if ~ismember(i,RemoveSet)
            for j = setdiff(i+1:length(TrainSet),RemoveSet)
                if dis(i,j) < sigma_radius
                    RemoveSet = [RemoveSet,j];
                end
            end
        end
    end
    Choose(RemoveSet) = false;
end