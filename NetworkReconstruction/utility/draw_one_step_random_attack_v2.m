function [A_obs, index_array, victim_likelihood_current, success] = draw_one_step_random_attack_v2(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array)

N = length(A_obs(1,:));
%     degree_array = zeros(1, N);
idx2ifr = sequence_missing(idx);
A_max = A_obs;
trials = 5;
trial = 1;
victim_likelihood_max = -inf;

while trial < trials
    A_temp = A_obs;
    index_array_temp = index_array;
    %Step 1: Draw sample from P(z_t-1| theta)
    for i = 1:N
        if i == idx2ifr
            continue
        end
        if index_array_temp(i) == 0
            continue
        end
        if A_temp(i, idx2ifr) == 1
            continue
        end
  
        dice = rand(1);
        if dice < P_k(psi(i), psi(idx2ifr))
            A_temp(i, idx2ifr) = 1;
            A_temp(idx2ifr, i) = 1;
        end 
    end
    index_array_temp(idx2ifr) = 1;
    %Calculate P(G_{t}|G_{t-1}, A)
    
    %3.Calculate the attack likelihood
    victim_likelihood = attack_likelihood_calc(A_temp, index_array_temp, idx2ifr, alpha);
    
    if victim_likelihood > victim_likelihood_max
        victim_likelihood_max = victim_likelihood;
        A_max = A_temp;
    end
    %index_array(idx2ifr) = 1;

    %victim_likelihood_current = 0;
    success = 1;
    trial = trial + 1;
end
victim_likelihood_current = victim_likelihood_max;
A_obs = A_max;
index_array(idx2ifr) = 1;