%This is a version that considers only the missing node
function [TP_sequence, FP_sequence, threshold, AUC, OPTROCPT] = AUC_calc_opt(A_obs_true, P_est_k, psi, sequence)
args = nargin;
true_labels = [];
pred_labels = [];

if args >3
    missing_hash = containers.Map('KeyType','int32','ValueType','int32');
    for i = 1:length(sequence)
        missing_hash(sequence(i)) = 1;
    end
    %Check on the adjacency matrix
    for i = 1: length(sequence)
        for j = 1: length(P_est_k(1,:))
            if sequence(i) == j
                continue
            end
            %         if isKey(missing_hash, j)
            %             continue
            %         end
            true_labels = [true_labels;A_obs_true(sequence(i), j)];
            pred_labels = [pred_labels;P_est_k(psi(sequence(i)), psi(j))];
        end
    end

    
else
    for i = 1:length(A_obs_true(1,:))-1
        for j = i+1:length(A_obs_true(1,:))
            true_labels = [true_labels;A_obs_true(i, j)];
            pred_labels = [pred_labels;P_est_k(psi(i), psi(j))];
            
        end
    end
end

[FP_sequence,TP_sequence,threshold,AUC, OPTROCPT] = perfcurve(true_labels,pred_labels,1);


