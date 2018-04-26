

function [A_obs, index_array, victim_likelihood, flag] = rejection_sampling_single_v2(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array)

if idx == length(sequence_missing) + 1
    victim_likelihood  = 0;
    flag = 1;
else
    trial = 0;
    %Backtracking when:
    %1. generate 1000 unsucessful samples or
    %2. can not generate a single sample after 1000 trials
    while trial < 10
        fprintf('At %d, try %d times \n', idx, trial);
        [A_obs, index_array, victim_likelihood_current, success] = draw_one_step_v4(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array);
        if success == 0
           flag = 0;
           victim_likelihood = 0;
           return
            
        else
            [A_obs, index_array, victim_likelihood_rest, success] = rejection_sampling_single_v2(A_obs, P_k, psi, idx + 1, sequence_missing, alpha, index_array);
            if success == 1
                victim_likelihood = victim_likelihood_current + victim_likelihood_rest;
                flag = 1;
                return
            else
                trial  = trial + 1;
                %Clean the adjacency matrix A_obs
                num_of_nodes = length(A_obs(1,:));
                for kk = 1:num_of_nodes
                    if A_obs(sequence_missing(idx),kk) == 1
                        A_obs(sequence_missing(idx),kk) = 0;
                        A_obs(kk, sequence_missing(idx)) = 0;
                    end
                end
                index_array(sequence_missing(idx)) = 0;
                continue
                
            end
        end     
    end
    victim_likelihood = 0;
    flag = 0;
end
    









