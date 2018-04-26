%Third method to draw a sample for one time step. 
%A MCMC method to take sample for one step


function [A_obs, index_array, victim_likelihood_current, success] = draw_one_step_v3(A_obs, P_k, psi, idx, sequence_missing, alpha, index_array)
trial = 0;
Burnin = 10;
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

if idx == 1 % the first node to infer
   ref_degree = avg_degree;    
else
   ref_degree = sum(A_obs(sequence_missing(idx-1), :)); 
end

prev_likelihood = 0;
prev_mu = ref_degree;
while 1
   accept = 0;
   %1. Draw a degree sample
    degree_array = zeros(1, N);
    %degree_sample = floor(ref_degree + rand(1)*avg_degree);
    degree_sample = floor(normrnd(prev_mu, sigma));
    counter = 0;
    temp_ll = 0;
%     while counter < degree_sample
%         if (I(kk) == idx2ifr || index_array(I(kk)) == 0)
%             kk = kk + 1;
%             continue
%         end
%         A_obs(idx2ifr, I(kk)) = 1;
%         A_obs(I(kk), idx2ifr) = 1; 
%         counter  = counter + 1;  
%         kk = kk + 1;
%     end
    %2.Update the adjacency matrix and calculate the cnnct ll
    for kk = 1:length(I)
        if I(kk) == idx2ifr
            continue
        end
        if counter < degree_sample
            temp_ll = temp_ll + log(P_k(psi(I(kk)), psi(idx2ifr)));
            A_obs(idx2ifr, I(kk)) = 1;
            A_obs(I(kk), idx2ifr) = 1; 
        else
            temp_ll = temp_ll + log(1-P_k(psi(I(kk)), psi(idx2ifr)));
        end
        counter = counter + 1;
    end
    
    %3.Calculate the attack likelihood
    norm_factor = 0;
    for i = 1:N
        degree_array(i) = sum(A_obs(i,:));
        norm_factor = norm_factor + degree_array(i)^(alpha);
    end
    victim_likelihood = degree_array(idx2ifr)^(alpha)/norm_factor;
    
    cur_likelihood = log(victim_likelihood) + temp_ll;
    
    
    
    %3.Accept or reject based on the previously drawn sample
    if trial == 0 % accept always the first sample
        prev_likelihood = cur_likelihood;
        trial = trial + 1;
        accept = 1;
        
    else % accept with prob A = min(1, p(cur)/p(prev))
        diff = cur_likelihood - prev_likelihood;
        if diff >= 0 % accept the sample
            prev_likelihood = cur_likelihood;
            accept = 1;
            prev_mu = sum(A_obs(idx2ifr,:));
        else
            accept_prob = exp(diff);
            dice = rand(1);
            if dice < accept_prob
                prev_likelihood = cur_likelihood;
                accept = 1;
                prev_mu = sum(A_obs(idx2ifr,:));
            end
        end
    end
    if trial >= Burnin && accept == 1
        victim_likelihood_current = log(victim_likelihood);
        success = 1;
        index_array(idx2ifr) = 1;
        degree_array(idx2ifr)
        return
    end
    %Clean the adjacency matrix for next sample
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
  
end

end