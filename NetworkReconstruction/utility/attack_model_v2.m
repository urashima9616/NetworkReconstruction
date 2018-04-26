function [G_new, sequence, index_array, max_sequence] = attack_model_v2(G, alpha, time_steps)
A = adjacency(G);
[N, m] = size(A);
index_array = ones(1,N);
degree_array = zeros(1, N);
sequence = zeros(1, time_steps);

for i = 1:N
    degree_array(i) = sum(A(i,:));
end


norm_factor = norm_factor_calc(degree_array, alpha);
%degree_array_augemented = degree_array.^(alpha);
%norm_factor = sum(degree_array_augemented);
max_sequence = zeros(2,time_steps);
for t = 1:time_steps
    %Generate the partition on [0,1]
    start = 0;
    partitions = zeros(1, N);
    for i = 1:N
        if index_array(i) == 0 || degree_array(i) == 0
            if i == 1
                partitions(i) = 0;
            else
                partitions(i) = partitions(i-1);
            end
        else
            delta = (degree_array(i)^(alpha))/norm_factor;
            partitions(i) = start + delta;
        end
        start = partitions(i); 
    end
    %Draw nodes from the current graph
    max_sequence(1,t) = max(degree_array);
    dice = rand(1);
    while dice == 0
        dice = rand(1);
    end
        
    for i = 1:N
        
        if dice <= partitions(i)
            victim_id = i;
            break
        end
    end
    max_sequence(2,t) = degree_array(victim_id);
    %Update sequence, index, adj matrix, degree and normalization factor
    sequence(t) = i;
    index_array(i) = 0;
    %norm_factor = norm_factor - degree_array(victim_id)^(alpha);
    degree_array(victim_id) = 0;
    for j = 1:N
        if A(victim_id, j) == 1
            A(victim_id, j) = 0;
            A(j, victim_id) = 0;
           %norm_factor = norm_factor - degree_array(j)^(alpha);
            degree_array(j) = degree_array(j) - 1;
           %norm_factor = norm_factor + degree_array(j)^(alpha);
        end
        
    end
    norm_factor = norm_factor_calc(degree_array, alpha);
    
end
%Generate new graph
G_new = graph(A);