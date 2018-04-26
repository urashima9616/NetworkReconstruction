%This is a generator for multi-fractal generative network
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
function [G, cor_index, m_opt, k_opt, P] = MFNG_gen(num_of_nodes, degree_avg, forced_connct)
N = num_of_nodes;
k_upper = 4;
m_upper = 10;

m_candidates = zeros(k_upper, m_upper);
for k = 1:k_upper
    for m = 1:m_upper
        if (m^(2*k)) <= N/degree_avg
            m_candidates(k,m) = N*m^(-2*k);
        end
    end                 
end

m_candidates = abs(degree_avg - m_candidates);
min_diff = degree_avg;
k_opt = 1;
m_opt = 1;
for i=1:k_upper
    for j = 1:m_upper
        if m_candidates(i,j) < min_diff
            k_opt = i;
            m_opt = j;
            min_diff = m_candidates(i,j);
        end
    end
end
            
actual_degree_avg = N*m_opt^(-2*k_opt);
num_of_parts = m_opt^k_opt;
expected_num_nodes = degree_avg*m_opt^(2*k_opt);
fprintf('The expected average degree is %d and the actual one is %f\n', degree_avg, actual_degree_avg);
fprintf('m = %d, k= %d\n', m_opt, k_opt);
fprintf('The expected number of nodes under this combination is %d\n', expected_num_nodes);
%randomly generate coordinates for N nodes
N = expected_num_nodes;%override the number of nodes, for debug use.
step = 1/num_of_parts;%equally sized boxes
cor_index = floor(rand(1, N)/step)+1;
cor_index(cor_index >= num_of_parts) = num_of_parts;


%Get the generative measure P

P = rand(m_opt);
P = triu(P) + triu(P)';
P = P./sum(sum(P));


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
    for j = i+1:N
        if (i == j)
            continue;
        elseif (adj_matrix(j,i) == 1  || adj_matrix(i,j) == 1)
            continue;
        else
            source_index = cor_index(i);
            dst_index = cor_index(j);
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