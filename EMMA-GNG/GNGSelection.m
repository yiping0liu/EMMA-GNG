function Parents = GNGSelection(finalPop,DiscretePop,Population,TrainSet,EdgeDis,net,NodeSolution,emptyIndex,Problem)      
        Pro = (Problem.FE-Problem.maxFE*.1)/(0.9*Problem.maxFE); 
        CellSet = setdiff(1:length(net.w),emptyIndex); 
        numberofnodes = length(CellSet);   
        Parent1s = [];
        Parent2s = [];
        
        NDis = pdist2(net.w(CellSet,:),net.w(CellSet,:));  
        NDis(logical(eye(numberofnodes))) = inf;  
        [sort_NDis,I] = sort(NDis);
        final_NDis = sort_NDis(1,:);
        Index_node = final_NDis > EdgeDis;     
        node1 = find(Index_node==1);        
       
        if length(DiscretePop) < 10
            firstNum = 0;
            seondNum = floor((1-Pro)*Problem.N);
        else
            firstNum = floor(0.5*(1-Pro)*Problem.N);
            seondNum = floor(0.5*(1-Pro)*Problem.N);
        end            
        ThirdNum = Problem.N - firstNum - seondNum;  
        
        [Parentstemp1,Parentstemp2] = expandSelection(finalPop,Population,net,NodeSolution,TrainSet,emptyIndex,seondNum);
        Parent1s = [Parent1s;Parentstemp1];
        Parent2s = [Parent2s;Parentstemp2];
        
        [Parentstemp1,Parentstemp2] = neighborSelection(finalPop,net,node1,NodeSolution,emptyIndex,ThirdNum);
        Parent1s = [Parent1s;Parentstemp1];
        Parent2s = [Parent2s;Parentstemp2];
       
        pDis = pdist2(DiscretePop.decs,DiscretePop.decs);
        pDis(logical(eye(length(DiscretePop)))) = inf;        
        [~,I] = sort(pDis);
		
        for i = 1:firstNum
            j = randperm(length(DiscretePop),1);
            Parent1 = DiscretePop(j).decs;
            Pool = I(randperm(10,1),j);
            Parent2 = DiscretePop(Pool).decs;
            Parent1s = [Parent1s;Parent1];
            Parent2s = [Parent2s;Parent2];
        end
        Parents = [Parent1s;Parent2s]; 
end