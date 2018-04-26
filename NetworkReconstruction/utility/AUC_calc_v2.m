%This is a version that considers only the missing node
function [TP_sequence, FP_sequence] = AUC_calc_v2(A_obs_true, P_est_k, psi, sequence)
args = nargin;
steps = 100;
max_p = max(max(P_est_k));
min_p = min(min(P_est_k));
step_len = (max_p-min_p)/steps;
threshold = min_p-step_len:step_len:max_p;

TP_sequence = zeros(length(threshold), 1);
FP_sequence = zeros(length(threshold), 1);
missing_hash = containers.Map('KeyType','int32','ValueType','int32');
for i = 1:length(sequence)
    missing_hash(sequence(i)) = 1;
end

for kk = 1:length(threshold)
    kk
    A_obs_sample = P_est_k <= threshold(kk);
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
            if isKey(missing_hash, j)
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
    end
    TP_sequence(kk,1) = sum(TPs)./sum(Ps);
    FP_sequence(kk,1) = sum(FPs)./sum(Ns);
end



