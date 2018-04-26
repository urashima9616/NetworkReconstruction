function [TP_sequence, FP_sequence] = AUC_calc(A_obs_true, P_est_k, psi, sequence)
args = nargin;
steps = 100;
max_p = max(max(P_est_k));
step_len = (max_p)/steps;
threshold = 0:step_len:max_p;

TP_sequence = zeros(length(threshold), 1);
FP_sequence = zeros(length(threshold), 1);
if args < 4
    %Calcuate all the AUC curve
    for kk = 1:length(threshold)
        kk
        A_obs_sample = P_est_k <= threshold(kk);
        TPs = 0;
        FPs = 0;
        Ps = 0;
        Ns = 0;
        for i=1:length(P_est_k(1,:))-1
            for j=i+1:length(P_est_k(1,:))
                if i == j
                    continue
                end
                
                if A_obs_true(i, j) == 1
                    Ps = Ps + 1;
                else
                    Ns = Ns + 1;
                end
                
                if A_obs_true(i,j) == 1 && A_obs_sample(psi(i), psi(j)) == 1
                    TPs = TPs + 1;
                    continue
                end
                
                if A_obs_true(i,j) == 0 && A_obs_sample(psi(i), psi(j)) == 1
                    FPs = FPs + 1;
                    continue
                end
            end
        end
        TP_sequence(kk,1) = TPs./Ps;
        FP_sequence(kk,1) = FPs./Ns;
    end
else
    
    for kk = 1:length(threshold)
        kk
        A_obs_sample = P_est_k <= threshold(kk);
        visited = containers.Map('KeyType','int32','ValueType','int32');
        TPs = zeros(length(sequence),1);
        FPs = zeros(length(sequence),1);
        Ps = zeros(length(sequence),1);
        Ns = zeros(length(sequence),1);
        %Iterate on the true adjacency_matrix
        for i = 1: length(sequence)
            for j = 1: length(P_est_k(1,:))
                if sequence(i) == j
                    continue
                end
                if isKey(visited, sequence(i))
                    continue
                end
                %Get the trues and falses for this node
                if A_obs_true(sequence(i), j) == 1
                    Ps(i) = Ps(i) + 1;
                else
                    Ns(i) = Ns(i) + 1;
                end
                
                %Check for TP or FP
                %Check for TP
                if A_obs_sample(psi(sequence(i)), psi(j)) == 1 && A_obs_true(sequence(i), j) == 1
                    TPs(i) = TPs(i) + 1;
                    continue;
                end
                
                if A_obs_sample(psi(sequence(i)), psi(j)) == 1 && A_obs_true(sequence(i), j) == 0
                    FPs(i) = FPs(i) + 1;
                    continue;
                end
                
            end
            visited(sequence(i)) = 1;
        end
        TP_sequence(kk,1) = sum(TPs)./sum(Ps);
        FP_sequence(kk,1) = sum(FPs)./sum(Ns);
    end
    
end


