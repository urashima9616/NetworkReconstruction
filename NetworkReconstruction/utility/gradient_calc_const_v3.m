%This is a utility to calculate the gradient the LL(theta)
%

function gradient_const = gradient_calc_const_v3(P_est, k)
m = length(P_est(1,:));
gradient_const = zeros(size(P_est));
for i=1:k
    if i == 1
        P_k_est = P_est;
    else
        P_k_est = kron(P_k_est, P_est);
    end
end

for i = 1:length(P_k_est(1,:))
    for j = i:length(P_k_est(1,:))
        t_count = zeros(size(P_est));
        i_labels = calculate_labels(m, k, i);
        j_labels = calculate_labels(m, k, j);
        for kk = 1:length(i_labels)
            t_count(i_labels(kk), j_labels(kk)) =  t_count(i_labels(kk), j_labels(kk)) + 1;
        end
        
        %Update the gradient matrix with the contribution of current edge
        for l = 1:length(gradient_const(1,:))
            for n = l : length(gradient_const(1,:))
                gradient_const(l,n) = gradient_const(l,n) + ((P_k_est(i,j)/(1-P_k_est(i,j)))*(-t_count(l,n)))/P_est(l,n);
                gradient_const(n,l) = gradient_const(l,n); 
            end
        end  
    end
end




