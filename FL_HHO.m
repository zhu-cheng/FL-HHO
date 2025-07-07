function [RabbitEnergy, RabbitLocation, Curve, RunningTime] = SFLA_HHO(N, T, lb, ub, dim, fobj, X)
global ArmPath;
global SFLAHHO_ArmPath;
global HHO_ArmPath;
global WOA_ArmPath;
global ArmPath1;

GlobalUC = 0; 
PreRabbitEnergy = inf;
Curve = zeros(1, T);
SearchThreshold = 3; 

[E, J] = EJ(N, T);
nPopMemeplex = 5; 
nMemeplex = floor(N / nPopMemeplex);
I = reshape(1:N, nMemeplex, []);
fitness = inf(N, 1);
for i = 1:size(X, 1)
    fitness(i) = fobj(X(i, :));
end
[X, fitness] = SortPopulation(X, fitness);
RabbitLocation = X(1, :);
RabbitEnergy = fitness(1);

%%

for i = 1:nMemeplex
    MemeplexBestLocation(i, :) = X(I(i, 1), :); 
    MemeplexBestEnergy(i) = fitness(I(i, 1)); 
    
    MemeplexWorstLocation(i, :) = X(I(i, nPopMemeplex), :);  
    MemeplexWorstEnergy(i) = fitness(I(i, nPopMemeplex));  
end
SflaStage = 0;
t = 0;
HHOInvalidCounter = 0;
SFLAInvalidCounter = 0;
while t < T    
    for i = 1:nMemeplex       
        minMemeplexPop = min(X(I(i, :), :));
        maxMemeplexPop = max(X(I(i, :), :));
        
        SFLAMemeplexPop = X(I(i, :), :);
        SFLAMemeplexFitness = fitness(I(i, :));
        
        HHOMemeplexPop = SFLAMemeplexPop;
        HHOMemeplexFitness = SFLAMemeplexFitness;       
        SimpleSFLA = 0;        
        for j = 1:nPopMemeplex          
            EscapingEnergy =  E(I(i, j), t + 1);
            JumpStrength = J(I(i, j), t + 1);         
            ReferenceLocation = RabbitLocation;
            meanX = mean(X);            
              if j < nPopMemeplex   
                if abs(EscapingEnergy) >= 1                  
                    q = rand();                    
                    if GlobalUC >= SearchThreshold
                        rand_Hawk_index = floor(N * rand() + 1);
                        X_rand = X(rand_Hawk_index, :);
                    else
                        rand_Hawk_index = floor(nPopMemeplex * rand() + 1);
                        X_rand = X(I(i, rand_Hawk_index), :);
                    end                   
                    if q < 0.5                        
                        HHOMemeplexPop(j, :) = X_rand - rand() * abs(X_rand - 2 * rand() * HHOMemeplexPop(j, :));
                    elseif q >= 0.5
                        HHOMemeplexPop(j, :) = (ReferenceLocation - meanX) - rand() * ((ub - lb) * rand + lb);
                    end
                elseif abs(EscapingEnergy) < 1                 
                    r = rand();
                    if r >= 0.5 && abs(EscapingEnergy) < 0.5
                        HHOMemeplexPop(j, :) = ReferenceLocation - EscapingEnergy * abs(ReferenceLocation - HHOMemeplexPop(j, :));
                    end
                    if r >= 0.5 && abs(EscapingEnergy) >= 0.5  
                        HHOMemeplexPop(j, :) = (ReferenceLocation - HHOMemeplexPop(j, :)) - EscapingEnergy * abs(JumpStrength * ReferenceLocation - HHOMemeplexPop(j, :));
                    end                   
                    if r < 0.5 && abs(EscapingEnergy) >= 0.5 
                        X1 = ReferenceLocation - EscapingEnergy * abs(JumpStrength * ReferenceLocation - HHOMemeplexPop(j, :));
                        if fobj(X1) < fobj(HHOMemeplexPop(j, :))
                            HHOMemeplexPop(j, :) = X1;
                        else
                            X2 = ReferenceLocation - EscapingEnergy * abs(JumpStrength * ReferenceLocation - HHOMemeplexPop(j, :)) + rand(1, dim) .* Levy(dim);
                            if (fobj(X2) < fobj(HHOMemeplexPop(j, :)))
                                HHOMemeplexPop(j, :) = X2;
                            end
                        end
                    end
                    if r < 0.5 && abs(EscapingEnergy) < 0.5
                        X1 = ReferenceLocation - EscapingEnergy * abs(JumpStrength * ReferenceLocation - mean(X(I(i, :), :)));
                        if fobj(X1) < fobj(HHOMemeplexPop(j, :)) 
                            HHOMemeplexPop(j, :) = X1;
                        else 
                            X2 = ReferenceLocation - EscapingEnergy * abs(JumpStrength * ReferenceLocation - mean(X(I(i, :), :))) + rand(1, dim) .* Levy(dim);
                            if (fobj(X2) < fobj(HHOMemeplexPop(j, :))) 
                                HHOMemeplexPop(j, :) = X2;
                            end
                        end
                    end
                end

                FU = HHOMemeplexPop(j, :) > ub;
                FL = HHOMemeplexPop(j, :) < lb;
                HHOMemeplexPop(j, :) = (HHOMemeplexPop(j, :) .* (~(FU + FL))) + (ub - rand()*power(10, -10)) .* FU + (lb + rand()*power(10, -10)) .* FL;                
                HHOMemeplexFitness(j) = fobj(HHOMemeplexPop(j, :));
            end           
            if SFLAInvalidCounter < 2 * nMemeplex
               [SFLAMemeplexPop, SFLAMemeplexFitness] = F(SimpleSFLA, fla_params, HHOMemeplexPop, HHOMemeplexFitness);          
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
            MemeplexPop = [MemeplexPop; HHOMemeplexPop];
            MemeplexFitness = [MemeplexFitness; HHOMemeplexFitness];
            
            X(I(i, :), :) =  HHOMemeplexPop;
            fitness(I(i, :)) = HHOMemeplexFitness;
            
        end
        %¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü¡ü        
    end  
    if SFLAInvalidCounter >= 2 * nMemeplex        
        SFLAInvalidCounter = 0;
    end   
    [X, fitness, SortOrder] = SortPopulation(X, fitness);       
    PreRabbitEnergy = RabbitLocation;    
    if fitness(1) < RabbitEnergy
        PreRabbitEnergy = RabbitLocation;
        RabbitLocation = X(1, :);
        RabbitEnergy = fitness(1);
        fobj(RabbitLocation);
    end
    t = t + 1;    
    Curve(t) = RabbitEnergy;
end

end




