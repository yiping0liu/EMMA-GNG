function [net,NodeSolution,emptyIndexFinal,finalPop,DiscretePop,EdgeDis] = ImproveGNG(net,Archive,Population,UEA,params,PlotFlagLabel)

    if ~exist('PlotFlag', 'var')
        PlotFlag = PlotFlagLabel;
    end

    Data = UEA.decs;
    DataSize = size(Data,1);
   
    %% Parameters
    N = params.N/2;   
    L = params.L;
    epsilon_b = params.epsilon_b;
    epsilon_n = params.epsilon_n;
    alpha = params.alpha;
    delta = params.delta;
    T = params.T;
    maxIter = 1;
    
    t = net.t;
    w = net.w;
    E = net.E;
    C = net.C;
    nx = net.nx;
    
    hp = net.hp;
    maxHP = 2*DataSize;  
    
    %% Randamize Input Data
    [nData,~] = size(Data); 
    Signals = Data(randperm(nData), :);
    
    %% Loop
    for iter = 1:maxIter      
        for kk = 1:nData            
            % Select Input
            nx = nx + 1;
            x = Signals(kk,:);
            
            % Competion and Ranking
            d = pdist2(x, w);
            [~, SortOrder] = sort(d);
            s1 = SortOrder(1);
            s2 = SortOrder(2);

            % change HP
            hp = hp-1; 
            hp(s1) = maxHP;    
            hp(s2) = hp(s2)+1;    
            
            
            t(s1, :) = t(s1, :) + 1;
            t(:, s1) = t(:, s1) + 1;

            E(s1) = E(s1) + d(s1)^2;

            w(s1,:) = w(s1,:) + epsilon_b*(x-w(s1,:));
            Ns1 = find(C(s1,:)==1);
            for j=Ns1
                w(j,:) = w(j,:) + epsilon_n*(x-w(j,:));
            end

            C(s1,s2) = 1;
            C(s2,s1) = 1;
            t(s1,s2) = 0;
            t(s2,s1) = 0;


            C(t>T) = 0;
            
            
            DeadNodes = (hp<=0);
            C(DeadNodes, :) = [];
            C(:, DeadNodes) = [];
            t(DeadNodes, :) = [];
            t(:, DeadNodes) = [];
            w(DeadNodes, :) = [];
            E(DeadNodes) = [];
            hp(DeadNodes) = [];
            
            if mod(nx, L) == 0 && size(w,1) < N   
                [~, q] = max(E);
                [~, f] = max(C(:,q).*E);
                r = size(w,1) + 1;
                w(r,:) = (w(q,:) + w(f,:))/2;
                C(q,f) = 0;
                C(f,q) = 0;
                C(q,r) = 1;
                C(r,q) = 1;
                C(r,f) = 1;
                C(f,r) = 1;
                t(r,:) = 0;
                t(:, r) = 0;
                E(q) = alpha*E(q);
                E(f) = alpha*E(f);
                E(r) = E(q);
                hp(r) = maxHP;                
            end

            E = delta*E;
            
        end
    end
    
    nNeighbor = sum(C);
    AloneNodes = (nNeighbor==0);
    nodeindex = find(AloneNodes==1);  
    nodeofnumber = sum(AloneNodes);   
    
    connection = graph(C~=0);    
    
    temp =  table2array(connection.Edges);
    dis = [];
    for i = 1:size(temp,1)
        dis = [dis;pdist2(w(temp(i,1),:),w(temp(i,2),:))];
    end
    EdgeDis = mean(dis);
    

    wholePop = [Archive,Population];
    [CombineData,ia] = unique(wholePop.decs,'rows','stable');
    wholePop = wholePop(ia);
    if EdgeDis > 0.1
        EdgeDis = 0.1;
    end
    thresholdDec = EdgeDis;    
    DataDisDec = pdist2(w,wholePop.decs);    
    [MinDataDisDec,~] = min(DataDisDec,[],1);
    NearIndex = find(MinDataDisDec <= thresholdDec);
    tempData = wholePop(NearIndex).decs;
    tempPop = wholePop(NearIndex);
   
    FarIndex = MinDataDisDec > thresholdDec; 
    DiscretePop = wholePop(FarIndex);
  
    NodeSolution = cell(1,size(w,1));
    NSDis = pdist2(w,tempData);  
    [~,NDIndex] = min(NSDis,[],1);
    for i = 1:length(w)
        NodeSolution{i} = find(NDIndex == i);
    end
    
    emptyIndex = [];
    for i=1:length(w)
        if isempty(NodeSolution{i})
            emptyIndex = [emptyIndex,i];
        end
    end
    
    Set = [];
    Dis = pdist2(w,Data);  
    [~,DIndex] = min(Dis,[],1);
    for i = emptyIndex
        tempIndex = find(DIndex == i)';
        Set = [Set;tempIndex];
    end
    finalPop = [tempPop,UEA(Set)];
    FDis = pdist2(w,finalPop.decs); 
    [~,FIndex] = min(FDis,[],1);
    for i = emptyIndex
        NodeSolution{i} = find(FIndex == i);
    end

    emptyIndexFinal = [];  
    for i=1:length(w)
        if isempty(NodeSolution{i})
            emptyIndexFinal = [emptyIndexFinal,i];
        end
    end

    if PlotFlag
        PlotResults(finalPop.decs, w, C);
        pause(0.01);
    end

    net.w = w;
    net.E = E;
    net.C = C;
    net.t = t;
    net.nx = nx;
    net.hp = hp;   
end