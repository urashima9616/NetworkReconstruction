function victim_likelihood = attack_likelihood_calc(A_obs, index_array, victim_idx, alpha)

%3.Calculate the attack likelihood
    N = length(A_obs(1,:));
    norm_factor = 0;
    degree_array = zeros(1,N);
    for i = 1:N
        if index_array(i) == 0
            continue
        end
        degree_array(i) = sum(A_obs(i,:));
        if degree_array(i) == 0
            continue
        end
        norm_factor = norm_factor + degree_array(i)^(alpha);
    end
    
    victim_likelihood =  (alpha)*log(degree_array(victim_idx)) - log(norm_factor);
   
end