function gradient = num_grad_calc(A_obs, P_est, P_k_est, psi, m, k)

delta = 0.001;
LL_original = LL_calc(A_obs, P_k_est, psi);

gradient = zeros(m);
for i = 1:m
    for j = 1:m
        P_est_sig = P_est;
        P_est_sig(i,j) = P_est_sig(i,j) + delta;
        for tt=1:k
            if tt == 1
                P_k_est_sig = P_est_sig;
            else
                P_k_est_sig = kron(P_k_est_sig, P_est_sig);
            end
        end
        LL_delta = LL_calc(A_obs, P_k_est_sig, psi);
        gradient(i,j) = (LL_delta - LL_original)./delta;
    end
end