function [Parent1s,Parent2s] = neighborSelection(finalPop,net,node,NodeSolution,emptyIndex,N)
    TrainPopulation = finalPop;
    Parent1s = [];
    Parent2s = [];
    CellSet = setdiff(1:length(net.w),emptyIndex);
    number = 1;  
    if isempty(node)
        cycles = 1:length(CellSet);
    else
        cycles = node;
    end
    remainNode = setdiff(1:length(CellSet),cycles);
    randIndex = randperm(length(cycles));
    cycles = cycles(randIndex);  
    Connection = distances(graph(net.C));
    for i = 1:ceil(N/length([cycles,cycles,remainNode]))   
        for j = [cycles,cycles,remainNode]
            Parent1_Pool = [NodeSolution{CellSet(j)}];     
            index = Parent1_Pool(randperm(end,1));
            Parent1 = TrainPopulation(index).decs;
            Parent1s = [Parent1s;Parent1];

            [~,conIndex] = find(Connection(CellSet(j),:) == 1); 
            Pool = [NodeSolution{[CellSet(j),conIndex]}]; 
            Pool = setdiff(Pool,index);
            if length(Pool) < 1
                remainPop = setdiff(1:length(TrainPopulation),index);  
                Dis = pdist2(net.w(CellSet(j),:),TrainPopulation(remainPop).decs);     
                    [~,I] = min(Dis,[],2);
                    Parent2 = TrainPopulation(remainPop(I)).decs;
            else
                P2 = Pool(randperm(end,1));
                Parent2 = TrainPopulation(P2).decs;
            end
            Parent2s = [Parent2s;Parent2];
            number = number +1;
            if number > N
                break;
            end
        end
    end   
end
