function [Population,fDKN] = ArchiveUpdate(Population,N,K)
% Archive Update in CPDEA

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Yen, G.G., Nojima, Y. and Masuyama, N., 2020. 
% Handling imbalance between convergence and diversity in the decision 
% space in evolutionary multimodal multiobjective optimization. IEEE 
% Transactions on Evolutionary Computation, 24(3), pp.551-565.
%------------------------------- Copyright --------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    FrontNo = NDSort(Population.objs,N);   
    Population = Population(FrontNo ==1);
    
     %% Remove duplicate solutions
    [~,ia] = unique(Population.decs,'rows');
    Population = Population(ia);
    
    [Choose,fDKN] = DoubleNearestSelection(Population.objs,Population.decs,N,K);
    
    Population = Population(Choose);
    fDKN = fDKN(Choose);
    
   
end
function [Choose,fDN] = DoubleNearestSelection(PopObj,PopDec,N,K)
% Select solutions based on double k-nearest neighbor method

    Np = size(PopObj,1);
    Choose = true(1,Np);
    
    if Np <= K
        fDN = zeros(1,Np);
        return;
    end   
    
    d_obj = pdist2(PopObj,PopObj,'euclidean');
    d_dec = pdist2(PopDec,PopDec,'euclidean');
%     d_obj(logical(eye(Np))) = inf;
%     d_dec(logical(eye(Np))) = inf;   % k=2 for better results
    
    sdo = sort(d_obj);
    sdd = sort(d_dec);    
    dn_obj = sum(sdo(1:K,:));      
    dn_dec = sum(sdd(1:K,:));
    avg_dn_obj = mean(dn_obj);
    avg_dn_dec = mean(dn_dec);
    if avg_dn_obj == 0
        avg_dn_obj = inf;
    end
    if avg_dn_dec == 0
        avg_dn_dec = inf;
    end
    fDN = 1./(1+dn_obj./avg_dn_obj+dn_dec./avg_dn_dec);
    
                
    while sum(Choose) > N
        [~,Del] = max(fDN);    %只找到第一个最大值
        Choose(Del) = false;
        d_obj(Del,:) = inf;
        d_obj(:,Del) = inf;
        d_dec(Del,:) = inf;
        d_dec(:,Del) = inf;
        
        sdo = sort(d_obj);
        sdd = sort(d_dec);    
        dn_obj = sum(sdo(1:K,:));
        dn_dec = sum(sdd(1:K,:));
        avg_dn_obj = mean(dn_obj(Choose));   % 更新一下只计算剩下的个体的均值
        avg_dn_dec = mean(dn_dec(Choose));
        if avg_dn_obj == 0
            avg_dn_obj = inf;
        end
        if avg_dn_dec == 0
            avg_dn_dec = inf;
        end
        fDN = 1./(1+dn_obj./avg_dn_obj+dn_dec./avg_dn_dec);
               
        fDN(~Choose) = -inf;
    end         
end

