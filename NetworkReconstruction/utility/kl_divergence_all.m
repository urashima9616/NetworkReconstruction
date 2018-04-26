function kl  = kl_divergence_all(P_est_k, P_k)

kl = 0;
for i = 1:length(P_est_k(1,:))-1
    for j = i+1:length(P_k(1,:))
        kl = kl -(P_k(i,j)* log(P_est_k(i,j)/P_k(i,j)));  
    end
end
