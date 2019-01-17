%This is a helper function to calculate the k-dimension labels given a
%index in P^k matrix

function [labels] = calculate_labels(m, k, index)
res = index;
labels = zeros(1, k);
len = k;
while k >0
    weight = m^(k-1);
    residual = mod(res, weight);
    result = floor(res/weight);
    
    if residual == 0
        labels(len-k+1) = result;
        res = weight;
    else
        labels(len-k+1) = result  + 1;
        res = res - weight*result;
    end
    k = k -1;
end
