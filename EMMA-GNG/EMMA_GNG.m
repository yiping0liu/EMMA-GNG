classdef EMMA_GNG < ALGORITHM
% <multi/many> <real> <multimodal>
% Evolutionary multimodal multiobjective optimization guided by growing 
% neural gas
% K --- 3 --- parameter for k-nearest
% eta --- 2 --- parameter for local convergence

%------------------------------- Reference --------------------------------
% Liu Y, Zhang L, Zeng X, et al. Evolutionary multimodal multiobjective 
% optimization guided by growing neural gas[J]. Swarm and Evolutionary 
% Computation, 2024, 86: 101500.
%------------------------------- Copyright --------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
   
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [K,eta] = Algorithm.ParameterSet(3,2);
            alpha = 0.1;
            %% GNG Parameter setting
            params.N = Problem.N;     % Max number of nodes
            params.MaxIt = 50;  % Maximum Generation
            params.L = 50;     % Cycle for topology reconstruction 
            params.epsilon_b = 0.2;     % Learning coefficient
            params.epsilon_n = 0.006;   % Learning coefficient of neighbor
            params.alpha = 0.5;     % Maximum cumulative error and its Maximum cumulative error neighbor reduction constant
            params.delta = 0.995;  % Error reduction coefficient
            params.T = 100;     %the value of amax
            
            %% Initialization
            Population = Problem.Initialization(Problem.N); % randomly generate N solutions    
            Archive = ArchiveUpdate(Population,Problem.N,K); % update archive     
            D_Dec = pdist2(Population.decs,Population.decs,'euclidean'); % distance between pair of solutions in the decision space
            DominationX = zeros(Problem.N); % Pareto domination relationship between pair of solutions
            for i = 1:Problem.N
                L1 = Population(i).objs < Population.objs;
                L2 = Population(i).objs > Population.objs;
                K1 = all(L1|(~L2),2);
                K2 = all(L2|(~L1),2);
                DominationX(i,K1) = 0;
                DominationX(K1,i) = 1;
                DominationX(i,K2) = 1;
                DominationX(K2,i) = 0;
            end
            [Population,fCPD,~] = CPDSelection(Population,K,eta,D_Dec,DominationX,Problem); % Environmental Selection

            AddGNGLabel = false;
            InitialGNGLabel = false;
            TrainSet = Archive;    
            
            %% Optimization
            while Algorithm.NotTerminated(Archive)             
                %% Generation Offsprings
                if  AddGNGLabel   % second reproduction operator
                    Parents = GNGSelection(finalPop,DiscretePop,Population,TrainSet,EdgeDis,net,NodeSolution,emptyIndex,Problem);
                    Offspring = OperatorGAhalf(Parents);
                    Offspring = SOLUTION(Offspring); 
                else % first reproduction operator               
                    MatingPool = TournamentSelection(2,Problem.N,fCPD);
                    Offspring = OperatorGA(Population(MatingPool));
                end

                %% Update D_Dec and DominationX  
                CombinePopulation = [Population,Offspring];
                D_Dec = pdist2(CombinePopulation.decs,CombinePopulation.decs,'euclidean'); % distance between pair of solutions in the decision space
                size = length(CombinePopulation);
                DominationX = zeros(size); % Pareto domination relationship between pair of solutions;              
                for i=1:size
                   L1 = CombinePopulation(i).objs < CombinePopulation.objs;
                   L2 = CombinePopulation(i).objs > CombinePopulation.objs;
                   K1 = all(L1|(~L2),2);
                   K2 = all(L2|(~L1),2);
                   DominationX(i,K1) = 0;
                   DominationX(K1,i) = 1;  
                   DominationX(i,K2) = 1;
                   DominationX(K2,i) = 0;
                end  
                %% Environmental Selection
                [Population,fCPD,~] = CPDSelection([Population,Offspring],K,eta,D_Dec,DominationX,Problem);
                %% Update Archive
                Archive = ArchiveUpdate([Archive,Offspring],Problem.N,K); 
                TrainSet = TrainSetUpdate(TrainSet,Archive,Population,Offspring);
                if Problem.FE >=  Problem.maxFE*alpha
                    if ~InitialGNGLabel 
                        TrainSet = Clearing(TrainSet);    
                        [net,NodeSolution,emptyIndex,finalPop,DiscretePop,EdgeDis] = InitialGNG(Archive,Population,TrainSet,params,true);         
                        AddGNGLabel = true;
                        InitialGNGLabel = true;
                    else
                        TrainSet = Clearing(TrainSet);
                        [net,NodeSolution,emptyIndex,finalPop,DiscretePop,EdgeDis] = ImproveGNG(net,Archive,Population,TrainSet,params,true);
                    end 
                end
            end
        end
    end
end