%This is a likelihood evaluation function only based on the network model
%Input arguments:
%   G_obs_t: the graph at a time t
%   P_k:   the linking probability
%   mapped_index:   the corresponding index id of G_t-1
%   alpha : the parameter of the attack model
%   index_array: indicator of missing nodes at their hypothesized true index
%                in original graph.
function log_likelihood = LL_calc(A, P_k, psi)
[N, m] = size(A);
log_likelihood = 0;

for i= 1:N-1
    for j = i : N
        if A(i, j)  > 0
            log_likelihood = log_likelihood + log(P_k(psi(i), psi(j)));
        else
            log_likelihood = log_likelihood + log(1- P_k(psi(i), psi(j)));
        end
        
    end
end