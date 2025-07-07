% Harris's CV  hawk optimizer: In this algorithm, Harris' hawks try to catch the rabbit.

% T: maximum iterations
% N: populatoin size
% Curve: Convergence curve
% RunningTime: 运行时间
%round(rand)*2-1 生成1或-1

function [RabbitEnergy, RabbitLocation, Curve, RunningTime] = SFLA_HHO(N, T, lb, ub, dim, fobj, X)

tic;

global ArmPath;
global SFLAHHO_ArmPath;
global HHO_ArmPath;
global WOA_ArmPath;
global ArmPath1;

% if nargin < 8
%     App = [];
% end


%% initialize the location and Energy of the rabbit
% RabbitLocation = zeros(1, dim);
% RabbitEnergy = inf;

GlobalUC = 0;  %寻优时，全局最优值不变的次数

%PreRabbitLocation = zeros(1, dim);
PreRabbitEnergy = inf;

Curve = zeros(1, T);

SearchThreshold = 3; %寻优门槛
%
%%
%↓↓↓↓↓↓↓↓↓↓计算HHO的E和J↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
E = zeros(N, T);
J = zeros(N, T);

RecoveryQty = floor(3 * rand(1, 1)) + 3;   %rq 体力恢复的次数(3, 5)
RecoveryPoint = ExtractePoints(1:T, RecoveryQty + 2);  %均匀选取RecoveryQty + 2个数据点，包括1和T。

E_Start = [2 * ones(size(X, 1), 1), zeros(size(X, 1), RecoveryQty)];

E_End = 0;

mi = 0.4;

for i = 1:size(X, 1)
    
    RecoveryPointCounter = 1; %RecoveryPoint计数器
    
    for t = 1:1:T
        if t == RecoveryPoint(RecoveryPointCounter) && t ~= T
            
            if RecoveryPointCounter <= RecoveryQty + 1
                
                if RecoveryPointCounter > 1 && RecoveryPointCounter <= RecoveryQty   %&& RecoveryPointCounter <= RecoveryQty + 1
                    aaa1 = E_Start(i, RecoveryPointCounter - 1);
                    aaa2 = -2 * (RecoveryPoint(RecoveryPointCounter + 1) * RecoveryPoint(RecoveryPointCounter+1)) / (T * T) + 2;  %抛物线
                    E_Start(i, RecoveryPointCounter) = (aaa1 - aaa2) * rand + aaa2;
                elseif RecoveryPointCounter == RecoveryQty + 1
                    E_Start(i, RecoveryPointCounter) = mean(E_Start(i, 1:RecoveryQty));
                end
                
                tt = t - RecoveryPoint(RecoveryPointCounter) + 1;
                
                E(i, t) = (E_Start(i, RecoveryPointCounter) - E_End) * exp(-(power(tt, 2) / power((mi * (RecoveryPoint(RecoveryPointCounter + 1) - RecoveryPoint(RecoveryPointCounter))), 2))) + E_End;
                
                RecoveryPointCounter = RecoveryPointCounter + 1;
            end
            
        else
            tt = t - RecoveryPoint(RecoveryPointCounter - 1) + 1;
            E(i, t) = (E_Start(i, RecoveryPointCounter - 1) - E_End) * exp(-(power(tt + 1, 2) / power((mi * (RecoveryPoint(RecoveryPointCounter) - RecoveryPoint(RecoveryPointCounter - 1))), 2))) + E_End;
        end
        
        if t < RecoveryPoint(RecoveryPointCounter - 1) + T / (RecoveryQty + 1) / 2
            aaa = rand * (E(i, RecoveryPoint(RecoveryPointCounter - 1)) - E(i, t));
        else
            aaa = rand * (E(i, RecoveryPoint(RecoveryPointCounter - 1)) / 2 - E(i, t));
        end
        
        E(i, t) = (round(rand) * 2 - 1) * E(i, t) + (round(rand) * 2 - 1) * aaa;
        if E(i, t) > E(i, RecoveryPoint(RecoveryPointCounter - 1))
            E(i, t) = E(i, RecoveryPoint(RecoveryPointCounter - 1));
        elseif  E(i, t) < -E(i, RecoveryPoint(RecoveryPointCounter - 1))
            E(i, t) = -E(i, RecoveryPoint(RecoveryPointCounter - 1));
        end
        
        E0 = 2 * rand() - 1; %-1<E0<1
        
        E(i, t) = E(i, t) * E0;
        
        if t <= RecoveryPoint(RecoveryQty + 1)
            J(i, t) = (-1 / RecoveryPoint(RecoveryQty + 1) * t + 2) * (1 - rand()); %递减
        else
            J(i, t) = 2 * (1 - rand()); %最后一搏
        end
        
    end
end
%↑↑↑↑↑↑↑↑↑↑计算HHO的E和J↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑




%%
%↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
nPopMemeplex = 5;                          % Memeplex Size
nMemeplex = floor(N / nPopMemeplex);        % Number of Memeplexes

I = reshape(1:N, nMemeplex, []);
%↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑


%%
%↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
MemeplexUC = zeros(nMemeplex, 1);  %每个Memeplex在寻优时，最优值不变的次数

MemeplexBestLocation = zeros(nMemeplex, dim);  %每个Memeplex，组内最优解的位置
MemeplexBestEnergy = inf(nMemeplex, 1);  %每个Memeplex，组内最优解的代价

MemeplexWorstLocation = zeros(nMemeplex, dim);  %每个Memeplex，组内最差解的位置
MemeplexWorstEnergy = -inf(nMemeplex, 1);  %每个Memeplex，组内最差解的代价

%SFLACounter = zeros(nMemeplex, 1);

% LocalOptimum = 0;
%↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑


%%
% FLA Parameters
fla_params.q = max(round(0.3 * nPopMemeplex), 2);   % Number of Parents
fla_params.alpha = 3;   % Number of Offsprings
fla_params.beta = 5;    % Maximum Number of Iterations
fla_params.sigma = 2;   % Step Size
fla_params.CostFunction = fobj;
% fla_params.VarMin = lb;
% fla_params.VarMax = ub;
fla_params.VarSize = [1 dim];
%=====================================================

%%
fitness = inf(N, 1);
for i = 1:size(X, 1)
    fitness(i) = fobj(X(i, :));
end

[X, fitness] = SortPopulation(X, fitness);

RabbitLocation = X(1, :);
RabbitEnergy = fitness(1);

%%

for i = 1:nMemeplex
    MemeplexBestLocation(i, :) = X(I(i, 1), :);  %每个Memeplex，组内最优解的位置
    MemeplexBestEnergy(i) = fitness(I(i, 1));  %每个Memeplex，组内最优解的代价
    
    MemeplexWorstLocation(i, :) = X(I(i, nPopMemeplex), :);  %每个Memeplex，组内最差解的位置
    MemeplexWorstEnergy(i) = fitness(I(i, nPopMemeplex));  %每个Memeplex，组内最差解的代价
end



%%
% B = zeros(nPopMemeplex, dim, nMemeplex);
% fitnessB = inf(nMemeplex, nPopMemeplex);

%SFLACounter = zeros(nMemeplex, 1) * fla_params.alpha;
SflaStage = 0; %SFLA休眠阶段


%%
t = 0;
%SFLACounter = 0;
HHOInvalidCounter = 0;
SFLAInvalidCounter = 0;
while t < T
    
    
    for i = 1:nMemeplex
        
        %%
        %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
        
        minMemeplexPop = min(X(I(i, :), :));
        maxMemeplexPop = max(X(I(i, :), :));
        
        SFLAMemeplexPop = X(I(i, :), :);
        SFLAMemeplexFitness = fitness(I(i, :));
        
        HHOMemeplexPop = SFLAMemeplexPop;
        HHOMemeplexFitness = SFLAMemeplexFitness;
        
        %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
        %%
        SimpleSFLA = 0;
        
        for j = 1:nPopMemeplex %SFLA把nPopMemeplex当params.beta用
            
            
            %%
            %%
            %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
            EscapingEnergy =  E(I(i, j), t + 1);
            JumpStrength = J(I(i, j), t + 1);
            %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
            
            ReferenceLocation = RabbitLocation;
            meanX = mean(X);
            
            %%
            %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
%             if HHOInvalidCounter >= HHOThreshold ...
%                     && SFLAInvalidCounter >= SFLAThreshold
%                 HHOInvalidCounter = 0;
%                 SFLAInvalidCounter = 0;
%             end
            %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
            
            
            %if HHOInvalidCounter < HHOThreshold
           % if HHOInvalidCounter < nMemeplex ...
                  %  && j < nPopMemeplex
              if j < nPopMemeplex   
                if abs(EscapingEnergy) >= 1
                    %% Exploration:
                    q = rand();
                    
                    %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
                    if GlobalUC >= SearchThreshold
                        rand_Hawk_index = floor(N * rand() + 1);
                        X_rand = X(rand_Hawk_index, :);
                    else
                        rand_Hawk_index = floor(nPopMemeplex * rand() + 1);
                        X_rand = X(I(i, rand_Hawk_index), :);
                    end
                    %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
                    
                    if q < 0.5
                        %if q < 0.33
                        HHOMemeplexPop(j, :) = X_rand - rand() * abs(X_rand - 2 * rand() * HHOMemeplexPop(j, :));
                    elseif q >= 0.5
                        % elseif q >= 0.33 && q<0.66
                        HHOMemeplexPop(j, :) = (ReferenceLocation - meanX) - rand() * ((ub - lb) * rand + lb);
                        %                 else
                        %                     HHOMemeplexPop(j, :) = unifrnd(minMemeplexPop, maxMemeplexPop);
                    end
                elseif abs(EscapingEnergy) < 1
                    %% Exploitation:
                    %% phase 1: surprise pounce (seven kills)
                    r = rand(); % probablity of each event
                    if r >= 0.5 && abs(EscapingEnergy) < 0.5 % Hard besiege
                        HHOMemeplexPop(j, :) = ReferenceLocation - EscapingEnergy * abs(ReferenceLocation - HHOMemeplexPop(j, :));
                    end
                    if r >= 0.5 && abs(EscapingEnergy) >= 0.5  % Soft besiege
                        HHOMemeplexPop(j, :) = (ReferenceLocation - HHOMemeplexPop(j, :)) - EscapingEnergy * abs(JumpStrength * ReferenceLocation - HHOMemeplexPop(j, :));
                    end
                    %% phase 2: performing team rapid dives (leapfrog movements)
                    if r < 0.5 && abs(EscapingEnergy) >= 0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                        X1 = ReferenceLocation - EscapingEnergy * abs(JumpStrength * ReferenceLocation - HHOMemeplexPop(j, :));
                        if fobj(X1) < fobj(HHOMemeplexPop(j, :)) % improved move?
                            HHOMemeplexPop(j, :) = X1;
                        else % hawks perform levy-based short rapid dives around the rabbit
                            X2 = ReferenceLocation - EscapingEnergy * abs(JumpStrength * ReferenceLocation - HHOMemeplexPop(j, :)) + rand(1, dim) .* Levy(dim);
                            if (fobj(X2) < fobj(HHOMemeplexPop(j, :))) % improved move?
                                HHOMemeplexPop(j, :) = X2;
                            end
                        end
                    end
                    if r < 0.5 && abs(EscapingEnergy) < 0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                        X1 = ReferenceLocation - EscapingEnergy * abs(JumpStrength * ReferenceLocation - mean(X(I(i, :), :)));
                        if fobj(X1) < fobj(HHOMemeplexPop(j, :)) % improved move?
                            HHOMemeplexPop(j, :) = X1;
                        else % Perform levy-based short rapid dives around the rabbit
                            X2 = ReferenceLocation - EscapingEnergy * abs(JumpStrength * ReferenceLocation - mean(X(I(i, :), :))) + rand(1, dim) .* Levy(dim);
                            if (fobj(X2) < fobj(HHOMemeplexPop(j, :))) % improved move?
                                HHOMemeplexPop(j, :) = X2;
                            end
                        end
                    end
                end
                
                %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
                FU = HHOMemeplexPop(j, :) > ub;
                FL = HHOMemeplexPop(j, :) < lb;
               % HHOMemeplexPop(j, :) = (HHOMemeplexPop(j, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
                HHOMemeplexPop(j, :) = (HHOMemeplexPop(j, :) .* (~(FU + FL))) + (ub - rand()*power(10, -10)) .* FU + (lb + rand()*power(10, -10)) .* FL;     
                
                HHOMemeplexFitness(j) = fobj(HHOMemeplexPop(j, :));
                %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
                
%                 if min(HHOMemeplexFitness) >= RabbitEnergy
%                     HHOInvalidCounter = HHOInvalidCounter + 1;
%                     %SFLAInvalidCounter = SFLAInvalidCounter + 1;
%                 else
%                    HHOInvalidCounter = 0;
%                     %SFLAInvalidCounter = 0;
%                 end
            end
            %%
            %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
            %SFLA
           
            if SFLAInvalidCounter < 2 * nMemeplex
             % if SFLAInvalidCounter < nMemeplex
                
                if SimpleSFLA == 0
                    SFLAOffspringsInvalidCounter = 0;
                    SFLAInnerSort = 1;
                    for k = 1:fla_params.alpha   %Number of Offsprings
                        %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
                        %SFLAFlag = 0;
                        
                        if SFLAInnerSort == 1
                            [SFLAMemeplexPop, SFLAMemeplexFitness, ~] = SortPopulation(SFLAMemeplexPop, SFLAMemeplexFitness);
                        else
                            SFLAInnerSort = 1;
                        end
                        
                        ImprovementStep2 = false;
                        Censorship = false;
                        %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
                        
                        Step = fla_params.sigma * rand(fla_params.VarSize) .* (SFLAMemeplexPop(1, :) - SFLAMemeplexPop(end, :));
                        NewSol1 = SFLAMemeplexPop(end, :) + Step;
                        
                        if IsInRange(NewSol1, lb, ub)
                            NewSol1_Cost = fobj(NewSol1);
                            
                            if NewSol1_Cost < SFLAMemeplexFitness(end)
                                SFLAMemeplexPop(end, :) = NewSol1;
                                SFLAMemeplexFitness(end) = NewSol1_Cost;
                                SFLAOffspringsInvalidCounter = 0;
                            else
                                ImprovementStep2 = true;
                            end
                        else
                            ImprovementStep2 = true;
                        end
                        
                        % SFLA Improvement Step 2
                        if ImprovementStep2
                            Step = fla_params.sigma * rand(fla_params.VarSize) .* (RabbitLocation - SFLAMemeplexPop(end, :));
                            NewSol2 = SFLAMemeplexPop(end, :) + Step;
                            
                            if IsInRange(NewSol2, lb, ub)
                                NewSol2_Cost = fobj(NewSol2);
                                
                                if NewSol2_Cost < SFLAMemeplexFitness(end)
                                    SFLAMemeplexPop(end, :) = NewSol2;
                                    SFLAMemeplexFitness(end) = NewSol2_Cost;
                                    SFLAOffspringsInvalidCounter = 0;
                                else
                                    Censorship = true;
                                end
                            else
                                Censorship = true;
                            end
                        end
                        % Censorship
                        if Censorship
                            NewSolCen = unifrnd(minMemeplexPop, maxMemeplexPop);
                            NewSolCen_Cost = fobj(NewSolCen);
                            %SFLAMemeplexPop(end, :) = unifrnd(minMemeplexPop, maxMemeplexPop);
                            %SFLAMemeplexFitness(end) = fobj(SFLAMemeplexPop(end, :));
                            
                            if NewSolCen_Cost >= SFLAMemeplexFitness(end)
                                SFLAOffspringsInvalidCounter = SFLAOffspringsInvalidCounter + 1;
                                SFLAInnerSort = 0; %下一次计算Offspring，不进行排序
                            else
                                SFLAOffspringsInvalidCounter = 0;
                            end
                            
                            SFLAMemeplexPop(end, :) = NewSolCen;
                            SFLAMemeplexFitness(end) = NewSolCen_Cost;
                        end
                    end
                    if SFLAOffspringsInvalidCounter >= 2 * fla_params.alpha
                        SimpleSFLA = 1;
                    end
                else %SimpleSFLA == 1
                    SFLAMemeplexPop(end, :) = unifrnd(minMemeplexPop, maxMemeplexPop);
                    SFLAMemeplexFitness(end) = fobj(SFLAMemeplexPop(end, :));
              %    [SFLAMemeplexPop, SFLAMemeplexFitness, DEBestEnergy, DEBestLocation] = DE(nPopMemeplex, dim, lb, ub, SFLAMemeplexPop, SFLAMemeplexFitness, fobj, RabbitEnergy, RabbitLocation)
                end
                %end
                %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
                
                
                
                
                
                
                %%
                
                %             if fitness(I(i, j)) < MemeplexBestEnergy(i)
                %                 %ForeMemeplexBestEnergy(i) = MemeplexBestEnergy(i);
                %                 MemeplexBestEnergy(i) = fitness(I(i, j));
                %                 MemeplexBestLocation(i, :) = HHOMemeplexPop(j, :);
                %             elseif fitness(I(i, j)) > MemeplexWorstEnergy(i)
                %                 MemeplexWorstEnergy(i) = fitness(I(i, j));
                %                 MemeplexWorstLocation(i, :) = HHOMemeplexPop(j, :);
                %             end
                
                
                if min(SFLAMemeplexFitness) >= RabbitEnergy
                    SFLAInvalidCounter = SFLAInvalidCounter + 1;
                else
                    SFLAInvalidCounter = 0;
                    %             MemeplexPop = [MemeplexPop; SFLAMemeplexPop];
                    %             SMemeplexFitness = [SMemeplexFitness; SFLAMemeplexFitness];
                    %
                    %             X(I(i, :), :) =  MemeplexPop;
                    %             fitness(I(i, :)) = SMemeplexFitness;
                end
            end
        end
        
        %%
        MemeplexPop = [];
        MemeplexFitness = [];
        if SFLAInvalidCounter < nMemeplex ...
                &&  HHOInvalidCounter < nMemeplex
            MemeplexPop = [HHOMemeplexPop; SFLAMemeplexPop];
            MemeplexFitness = [HHOMemeplexFitness; SFLAMemeplexFitness];
            [MemeplexPop, MemeplexFitness, ~] = SortPopulation(MemeplexPop, MemeplexFitness);
            
            X(I(i, :), :) =  MemeplexPop(1:nPopMemeplex, :);
            fitness(I(i, :)) = MemeplexFitness(1:nPopMemeplex);
            
        elseif SFLAInvalidCounter < nMemeplex
            X(I(i, :), :) =  SFLAMemeplexPop;
            fitness(I(i, :)) = SFLAMemeplexFitness;
        elseif HHOInvalidCounter < nMemeplex
            %             if min(HHOMemeplexFitness) >= RabbitEnergy
            %                 HHOInvalidCounter = HHOInvalidCounter + 1;
            %             else
            %                 HHOInvalidCounter = 0;
            MemeplexPop = [MemeplexPop; HHOMemeplexPop];
            MemeplexFitness = [MemeplexFitness; HHOMemeplexFitness];
            
            X(I(i, :), :) =  HHOMemeplexPop;
            fitness(I(i, :)) = HHOMemeplexFitness;
            
        end
        %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑        
    end
    
%     if SFLAInvalidCounter < nMemeplex
%          SFLAInvalidCounter = 0;
%     end
%     
%     if HHOInvalidCounter < nMemeplex
%          HHOInvalidCounter = 0;
%     end
%     
%     if HHOInvalidCounter >= nMemeplex ...        
%          && SFLAInvalidCounter >= 2 * nMemeplex
%         HHOInvalidCounter = 0;
%         SFLAInvalidCounter = 0;
%     end
%     
    if SFLAInvalidCounter >= 2 * nMemeplex        
        SFLAInvalidCounter = 0;
    end
    %%
    
    [X, fitness, SortOrder] = SortPopulation(X, fitness);
    
    
    PreRabbitEnergy = RabbitLocation;
    
    if fitness(1) < RabbitEnergy
        PreRabbitEnergy = RabbitLocation;
        RabbitLocation = X(1, :);
        RabbitEnergy = fitness(1);
        fobj(RabbitLocation);
    SFLAHHO_ArmPath = ArmPath;
    end
    
%     SFLAHHO_ArmPath(t+1, 1:dim) = ArmPath;
%     SFLAHHO_ArmPath(t+1, dim+1) = fitness(1);
    %%
    t = t + 1;
    
    Curve(t) = RabbitEnergy;
    
    
%     if (PreRabbitEnergy - RabbitEnergy) / PreRabbitEnergy <= 0.0000001
%         GlobalUC = GlobalUC + 1;  %寻优时，全局最优值不变的次数
%     end
    
%     if HHOInvalidCounter >= nMemeplex ...
%             && SFLAInvalidCounter >= 2 * nMemeplex
%         
%         
%         [X, fitness, DEBestEnergy, DEBestLocation] = DE(N, dim, lb, ub, X, fitness, fobj, RabbitEnergy, RabbitLocation);
%         
%         if DEBestEnergy < RabbitEnergy
%             PreRabbitEnergy = RabbitLocation;
%             RabbitLocation = DEBestLocation;
%             RabbitEnergy = DEBestEnergy;
%         end
%         HHOInvalidCounter = 0;
%         SFLAInvalidCounter = 0;
%         
%     end
end
RunningTime = toc;
end






%%
%Levy飞行
function o = Levy(d)
beta = 1.5;
sigma = (gamma(1 + beta) * sin(pi * beta / 2)/(gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
u = randn(1, d) * sigma;
v = randn(1, d);
s = u ./ abs(v) .^ (1/beta);
o = s;
end

%%
%种群排序
%function [pop, SortOrder] = SortPopulation(pop)
function [X, fitness, SortOrder] = SortPopulation(X, fitness)
% Get Costs
%Costs = [pop.Cost];

% Sort the Costs Vector
%[~, SortOrder]=sort(Costs);
[f, i] = sort(fitness);

fitness = f;
X = X(i, :);
SortOrder = i;
% Apply the Sort Order to Population
%pop = pop(SortOrder);
end


function [pop, fitness, DEBestEnergy, DEBestLocation] = DE(N, dim, lb, ub, pop, fitness, fobj, RabbitEnergy, RabbitLocation)
beta_min=0.2;   % Lower Bound of Scaling Factor
beta_max=0.8;   % Upper Bound of Scaling Factor
pCR=0.2;        % Crossover Probability
DEBestEnergy = inf;
DEBestLocation = zeros(1, dim);
for i = 1:N
    x = pop(i, :);
    
    A = randperm( N );%1~N的随机排列
    
    A(A == i) = [];  %删掉i
    
    a = A( 1 );
    b = A( 2 );
    c = A( 3 );
    
    % Mutation
    beta = unifrnd(beta_min, beta_max, 1, dim);
    y = pop(a, :) + beta .* (pop(b, :) - pop(c, :));
    y = max(y, lb);
    y = min(y, ub);
    
    % Crossover
    z = zeros( size( x ) );
    j0 = randi([1 numel( x )]);  %产生1个1~numel(x)的随机整数
    for j = 1:numel( x )
        if j == j0 || rand <= pCR
            z( j ) = y( j );
        else
            z( j ) = x( j );
        end
    end
    
    
    NewSol = z;
    NewSol_Cost = fobj(NewSol);
    
    if NewSol_Cost < fitness(i)
        pop(i, :) = NewSol;
        fitness(i) = NewSol_Cost;
        
        if NewSol_Cost < DEBestEnergy
            DEBestLocation = NewSol;
            DEBestEnergy  = NewSol_Cost;
        end
    end
    
end
end


