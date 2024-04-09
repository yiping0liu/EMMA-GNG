function [Parent1s,Parent2s] = expandSelection(finalPop,Population,net,NodeSolution,TrainSet,emptyIndex,N)
    TotalPopulation = [finalPop,Population,TrainSet];
    [~,ia] = unique(TotalPopulation.decs,'rows','stable'); 
    TotalPopulation = TotalPopulation(ia);
    Parent1s = [];
    Parent2s = [];
    value = pi/6;
    CellSet = setdiff(1:length(net.w),emptyIndex);
    number = 0;  
    Connection = distances(graph(net.C));
    numberofnode = length(CellSet);
    Dec = pdist2(net.w(CellSet,:),TotalPopulation.decs);
    [~,I] = sort(Dec,2);
    while number < N
        j = randperm(numberofnode,1);
        Parent1_Pool = [NodeSolution{CellSet(j)}];      
        index = Parent1_Pool(randperm(end,1));
        Parent1 = TotalPopulation(index).decs;
        [~,conIndex] = find(Connection(CellSet(j),:) == 1); 
        cosine = 1 - pdist2(net.w(conIndex,:)-net.w(CellSet(j),:),TotalPopulation.decs-net.w(CellSet(j),:),'cosine');
        angle = acos(cosine);
        [~,c] = find(angle < value);
        neigbor = c';
        if size(c,1) == 1
            neigbor = c;
        end
        Pool = setdiff(I(j,:),[neigbor,Parent1_Pool],'stable');
        if length(Pool) < 10  
            continue;   
        else
            Pool = Pool(1:10);
            P2 = Pool(randperm(end,1));
        end
        Parent2 = TotalPopulation(P2).decs;
        Parent1s = [Parent1s;Parent1];
        Parent2s = [Parent2s;Parent2];
        number = number +1;
        if number > N
            break;
        end
    end   
end
