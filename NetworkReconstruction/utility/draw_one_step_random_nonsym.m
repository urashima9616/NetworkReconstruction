function [A_obs, index_array, victim_likelihood_current, success] = draw_one_step_random_nonsym(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array)

    N = length(A_obs(1,:));
%     degree_array = zeros(1, N);
    idx2ifr = sequence_missing(idx);
    %Step 1: Draw sample from P(z_t-1| theta)
    for i = 1:N
        if i == idx2ifr
            continue
        end
        %if index_array(i) == 0
        %    continue
        %end
        if A_obs(i, idx2ifr) == 1
            continue
        end
  
        dice = rand(1);
        if psi(i)>psi(idx2ifr)
            prob = P_k(psi(idx2ifr), psi(i));
        else
            prob = P_k(psi(i), psi(idx2ifr));
        end
        if dice < prob
            A_obs(i, idx2ifr) = 1;
            A_obs(idx2ifr, i) = 1;
        end 
    end
    %Calculate P(G_{t}|G_{t-1}, A)
    index_array(idx2ifr) = 1;

    victim_likelihood_current = 0;
    success = 1;
end