function [Population,fCPD,LocalC] = CPDSelection(Population,K,eta,D_Dec,DominationX,Prob)
% CPD based environmental selection in CPDEA

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Yen, G.G., Nojima, Y. and Masuyama, N., 2020. 
% Handling imbalance between convergence and diversity in the decision 
% space in evolutionary multimodal multiobjective optimization. IEEE 
% Transactions on Evolutionary Computation, 24(3), pp.551-565.
%------------------------------- Copyright --------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%-------------------------------------------------------------------------
    [Population,fCPD,LocalC]  = CPD(Population,K,eta,D_Dec,DominationX,Prob);
    %% remove the worst solution
    number = size(Population,2);
    while number > Prob.N
        [~,Del] = max(fCPD);
        Population(Del) = []; 
        D_Dec(Del,:) = [];
        D_Dec(:,Del) = [];
        DominationX(Del,:) = [];
        DominationX(:,Del) = [];
        [Population,fCPD,LocalC]  = CPD(Population,K,eta,D_Dec,DominationX,Prob);
        number = number-1;
    end
end

function [Population,fCPD,LocalC]  = CPD(Population,K,eta,D_Dec,DominationX,Prob)

    %% Local Convergence Quality 
    N = length(Population);
    V = prod(Prob.upper - Prob.lower);
    R = eta.*(V./N).^(1./Prob.D);
    temp = exp(-D_Dec.^2./(2.*R.^2))./(R.*(2.*pi).^.5); % Normal_distribution
    LocalC = NaN(1,N);
    for i = 1:N
        LocalC(i) = DominationX(i,:)*temp(:,i); 
    end

    %% Transform D_Dec based on Local Convergence Quality
    D_t = D_Dec;
    for i = 1:N-1
        x = (LocalC(i) + LocalC(i+1:end)).*0.5; 
        y = 1./((x+1).^1);
        D_t(i,i+1:end) = D_t(i,i+1:end).*y; 
        D_t(i+1:end,i) = D_t(i,i+1:end);
    end
    D_t(logical(eye(N))) = inf;
    %% k-nearest neighbor
    D_s = sort(D_t);
    fCPD = 1./(1+sum(D_s(1:K,:)));
end
