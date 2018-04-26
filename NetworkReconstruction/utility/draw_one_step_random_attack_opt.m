function [A_obs, index_array, victim_likelihood_current, success, norm_factor_max, degree_array_max] = draw_one_step_random_attack_opt(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array, norm_factor, degree_array)

N = length(A_obs(1,:));
%     degree_array = zeros(1, N);
idx2ifr = sequence_missing(idx);
A_max = A_obs;
degree_array_max = degree_array;
norm_factor_max = norm_factor;
trials = 4;
trial = 1;
victim_likelihood_max = -inf;

while trial < trials
    A_temp = A_obs;
    index_array_temp = index_array;
    norm_factor_temp = norm_factor;
    degree_array_temp = degree_array;
    %Step 1: Draw sample from P(z_t-1| theta)
    while degree_array_temp(idx2ifr) ==0
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
                %update a_obs
                A_temp(i, idx2ifr) = 1;
                A_temp(idx2ifr, i) = 1;
                %update norm_factor
                prev = degree_array(i).^(alpha);
                after =  (degree_array(i)+1).^(alpha);
                norm_factor_temp = norm_factor_temp - prev  + after;
                degree_array_temp(i) = degree_array_temp(i) + 1;
                degree_array_temp(idx2ifr) = degree_array_temp(idx2ifr) + 1;
            end
        end
        %degree_array_temp(idx2ifr);
    end
    share = degree_array_temp(idx2ifr).^(alpha);
    norm_factor_temp = norm_factor_temp +  share;

    index_array_temp(idx2ifr) = 1;
    %Calculate P(G_{t}|G_{t-1}, A)
    
    %3.Calculate the attack likelihood
    if share == 0
        victim_likelihood = -inf;
    else 
        victim_likelihood = log(share) - log(norm_factor_temp);
    end
    %attack_likelihood_calc(A_temp, index_array_temp, idx2ifr, alpha);
    
    if victim_likelihood > victim_likelihood_max
        victim_likelihood_max = victim_likelihood;
        A_max = A_temp;
        norm_factor_max =  norm_factor_temp;
        degree_array_max = degree_array_temp;
    end
    %index_array(idx2ifr) = 1;

    %victim_likelihood_current = 0;
    success = 1;
    trial = trial + 1;
end
victim_likelihood_current = victim_likelihood_max;
A_obs = A_max;
index_array(idx2ifr) = 1;
degree_array_max(idx2ifr)