%This is a one-step rejection sampler to take a single sample from P(z_t-1|theta, A, G_obs_t)
%Input arguments:
%   G_obs_t: the graph at a time t
%   P_k:   the linking probability
%   psi:  mapping between the index of G to index of P_k
%   idx2ifr:  index in G to infer
%   alpha : the parameter of the attack model
%   index_array: indicator of a node being a missing node or existing node
%Returns:
%G_sample: a sample of G_{t-1}
%index_array: updated index_array 
%victim_likelihood: the probablity that the inferred z_{t-1} being under
%attack given G_{t-1} = G_{t} + z_{t-1}

function [G_sample, index_array, victim_likelihood] = rejection_sampling_v2(G_obs, P_k, psi, sequence_missing, alpha, index_array)

while 1
   
    A = adjacency(G_obs);
    [N, m] = size(A);
    degree_array = zeros(1, N);
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
            A(i, idx2ifr) = 1;
            A(idx2ifr, i) = 1;
        end 
    end
    %Calculate P(G_{t}|G_{t-1}, A)
    index_array(idx2ifr) = 1;
    norm_factor = 0;

    for i = 1:N
        degree_array(i) = sum(A(i,:));
        norm_factor = norm_factor + degree_array(i)^(alpha);
    end
    victim_likelihood = degree_array(idx2ifr)^(alpha)/norm_factor;
    
    %Step 2: Draw uniformly from [0, P(z_t|theta)]
    %L = l_eval_log(A, P_k, psi, index_array);
    dice = rand(1);
    if dice < victim_likelihood
        break
    end
    
    %Step 3: Accept or Reject the sample
    %if dice <= L * victim_likelihood
    %    break
    %end
    index_array(idx2ifr) = 0;

end
G_sample = graph(A);







