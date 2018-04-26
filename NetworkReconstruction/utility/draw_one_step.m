function [A_obs, index_array, victim_likelihood_current, success] = draw_one_step(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array)
trial = 0;

while trial < 20
  
    N = length(A_obs(1,:));
    degree_array = zeros(1, N);
    idx2ifr = sequence_missing(idx);
    %Step 1: Draw sample from P(z_t-1| theta)
    for i = 1:N
        if i == idx2ifr
            continue
        end
        if index_array(i) == 0
            continue
        end
        dice = rand(1);
        if dice < P_k(psi(i), psi(idx2ifr))
            A_obs(i, idx2ifr) = 1;
            A_obs(idx2ifr, i) = 1;
        end 
    end
    %Calculate P(G_{t}|G_{t-1}, A)
    index_array(idx2ifr) = 1;
    norm_factor = 0;

    for i = 1:N
        degree_array(i) = sum(A_obs(i,:));
        norm_factor = norm_factor + degree_array(i)^(alpha);
    end
    victim_likelihood = degree_array(idx2ifr)^(alpha)/norm_factor;
    
    %Step 2: Draw uniformly from [0, P(z_t|theta)]
    %L = l_eval_log(A, P_k, psi, index_array);
    dice = rand(1);
    if dice < victim_likelihood
        victim_likelihood_current = log(victim_likelihood);
        success = 1;
        return
    else
        for i =1:N
            if A_obs(idx2ifr, i) == 1
                A_obs(idx2ifr, i) = 0;
                A_obs(i, idx2ifr) = 0;
            end
        end
        trial = trial + 1;
        
    end

    index_array(idx2ifr) = 0;

end

victim_likelihood_current = 0;
success = 0;
index_array(idx2ifr) = 0;
end