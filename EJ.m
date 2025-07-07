function [E, J] = EJ(N, T)
E = zeros(N, T);
J = zeros(N, T);

RecoveryQty = floor(3 * rand(1, 1)) + 3;   
RecoveryPoint = ExtractePoints(1:T, RecoveryQty + 2); 

E_Start = [2 * ones(size(X, 1), 1), zeros(size(X, 1), RecoveryQty)];

E_End = 0;

mi = 0.4;

for i = 1:size(X, 1)
    
    RecoveryPointCounter = 1;
    
    for t = 1:1:T
        if t == RecoveryPoint(RecoveryPointCounter) && t ~= T
            
            if RecoveryPointCounter <= RecoveryQty + 1
                
                if RecoveryPointCounter > 1 && RecoveryPointCounter <= RecoveryQty   
                    aaa1 = E_Start(i, RecoveryPointCounter - 1);
                    aaa2 = -2 * (RecoveryPoint(RecoveryPointCounter + 1) * RecoveryPoint(RecoveryPointCounter+1)) / (T * T) + 2;
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
        
        E0 = 2 * rand() - 1;
        
        E(i, t) = E(i, t) * E0;
        
        if t <= RecoveryPoint(RecoveryQty + 1)
            J(i, t) = (-1 / RecoveryPoint(RecoveryQty + 1) * t + 2) * (1 - rand()); 
        else
            J(i, t) = 2 * (1 - rand()); 
        end
        
    end
end
