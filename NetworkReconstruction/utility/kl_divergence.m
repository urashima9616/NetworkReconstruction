function kl  = kl_divergence(P_est_k, P_k, sequence)

kl = 0;

visited = containers.Map('KeyType','int32','ValueType','int32');

for i = 1:length(sequence)
    for j = 1:length(P_k(1,:))
        if sequence(i) == j
            continue
        end
        tf = isKey(visited, j);
        if tf
            continue
        end

        kl = kl -(P_k(sequence(i),j)* log(P_est_k(sequence(i),j)/P_k(sequence(i),j)));  
    end
    visited(sequence(i)) = 1;
end
