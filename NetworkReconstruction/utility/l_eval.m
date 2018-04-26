%This is a likelihood evaluation function only based on the network model
%Input arguments:
%   G_obs_t: the graph at a time t
%   P_k:   the linking probability
%   mapped_index:   the corresponding index id of G_t-1
%   alpha : the parameter of the attack model
%   index_array: indicator of missing nodes at their hypothesized true index
%                in original graph.
function likelihood = l_eval(A, P_k, psi, index_array)
[N, m] = size(A);
likelihood = 1;

for i= 1:N
    if index_array(i) == 0 % a missing node
        continue % ignore it
    end
    for j = i+1 : N
        if index_array(j) == 0 % a missing node
            continue % skip it
        end
        if A(i, j) == 1
            likelihood = likelihood * P_k(psi(i), psi(j));
        else
            likelihood = likelihood *  (1- P_k(psi(i), psi(j)));
        end
        
    end
end