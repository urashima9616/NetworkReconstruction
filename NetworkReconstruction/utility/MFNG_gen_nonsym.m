%This is a generator for multi-fractal generative network(Kronecker)
%Author: Yuankun Xue
%--------------------------------------------------------
%--Input arguments: 
%  num_of_nodes : Number of nodes to be generated. 
%  degree_avg: Expected average degree
%  forced_cnnct: The generated network wil be enforced to be connected if
%              non-zero value is set This is achieved by two ways:
%1.Reject disconnected network (forced_cnnct == 1)
%2.Normalize the link probability in a row-wise way as to make sure a node
%  is alway connected to at least one node (forced_cnnct == 2).
%  link_prob: Generating probability measure
%  parts 
%
%
%--Update History:
%12/24/2017:  Originally created by Yuankun Xue
function [G] = MFNG_gen_nonsym(k, m, P)
N = m^(k);
m_opt = m;
k_opt = k;


for i=1:k_opt
    if i == 1
        P_k = P;
    else
        P_k = kron(P_k, P);
    end
end
    
%Get the sample network
adj_matrix = zeros(N, N);
for i = 1:N
    for j = 1:N
        if (i == j)
            continue;
        elseif (adj_matrix(j,i) == 1  || adj_matrix(i,j) == 1)
            continue;
        else
            source_index = i;
            dst_index = j;
            dice = rand(1);
            if dice < P_k(source_index, dst_index)
                adj_matrix(i,j) = 1;
                adj_matrix(j,i) = 1;
            end
        end
    end 
end

G = graph(adj_matrix);
end