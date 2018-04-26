function [A_obs, index_array, victim_likelihood_current, success] = draw_one_step_v2(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array)
trial = 0;
N = length(A_obs(1,:));
idx2ifr = sequence_missing(idx);
[B, I] = sort(P_k(psi(idx2ifr),:));
I = fliplr(I);
avg_degree = sum(P_k(idx2ifr, :));

if idx == 1 % the first node to infer
   ref_degree = sum(sum(A_obs))/N;    
else
   ref_degree = sum(A_obs(sequence_missing(idx-1), :)); 
end


while trial < 100
  
    
    degree_array = zeros(1, N);
    %Step 1: Draw sample from A_alpha first
    
    degree_sample = floor(ref_degree + rand(1)*(avg_degree));
    %Distribute the links in descending order
    kk = 1;
    counter = 0;
    while counter < degree_sample
        if (I(kk) == idx2ifr || index_array(I(kk)) == 0)
            kk = kk + 1;
            continue
        end
        A_obs(idx2ifr, I(kk)) = 1;
        A_obs(I(kk), idx2ifr) = 1; 
        counter  = counter + 1;  
        kk = kk + 1;
    end
    
    
    %Calc the likelihood
    norm_factor = 0;
    for i = 1:N
        degree_array(i) = sum(A_obs(i,:));
        norm_factor = norm_factor + degree_array(i)^(alpha);
    end
    victim_likelihood = degree_array(idx2ifr)^(alpha)/norm_factor;
    
    %Accept or reject
    dice = rand(1);
    if dice < victim_likelihood
        victim_likelihood_current = log(victim_likelihood);
        success = 1;
        index_array(idx2ifr) = 1;
        degree_array(idx2ifr)
        
        return
    else
        kk = 1;
        counter = 0;
        while counter < degree_sample
            if (I(kk) == idx2ifr || index_array(I(kk)) == 0)
                kk = kk + 1;
                continue
            end
            A_obs(idx2ifr, I(kk)) = 0;
            A_obs(I(kk), idx2ifr) = 0; 
            counter  = counter + 1;  
            kk = kk + 1;
        end
        trial = trial + 1;
        
    end
end

victim_likelihood_current = 0;
success = 0;
index_array(idx2ifr) = 0;
end