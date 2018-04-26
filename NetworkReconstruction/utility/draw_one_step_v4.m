%4-th method to draw a sample for one time step. 
%Sample a degree from normal(avg_dergee, sigma) and accept or reject based
%on attack model


function [A_obs, index_array, victim_likelihood_current, success] = draw_one_step_v4(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array)
trial = 0;
N = length(A_obs(1,:));
idx2ifr = sequence_missing(idx);
[B, I] = sort(P_k(psi(idx2ifr),:));
I = fliplr(I);
degree_array_pre = zeros(1, N);
avg_degree = sum(P_k(idx2ifr, :));

for i =1 : N
   degree_array_pre(i) = sum(A_obs(i,:)); 
end
sigma = sqrt(var(degree_array_pre));

% if idx == 1 % the first node to infer
%    ref_degree = avg_degree;    
% else
%    ref_degree = sum(A_obs(sequence_missing(idx-1), :)); 
% end

% prev_likelihood = 0;
% prev_mu = ref_degree;
while trial < 100
   %1. Draw a degree sample
    degree_array = zeros(1, N);
    %degree_sample = floor(ref_degree + rand(1)*avg_degree);
    degree_sample = floor(normrnd(avg_degree, sigma));
    counter = 0;
    kk = 1;
    %2.Update the adjacency matrix and calculate the cnnct ll
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
    
    %3.Calculate the attack likelihood
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