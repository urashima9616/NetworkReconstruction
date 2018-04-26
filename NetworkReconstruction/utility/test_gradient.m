P_est_sig = P_est;
delta = 0.0001;
P_est_sig(2,2) = P_est_sig(2,2) + delta;
%A_obs = Z_samples{1};
G_sample = graph(A_obs);
for i=1:k_opt
    if i == 1
        P_k_est_sig = P_est_sig;
    else
        P_k_est_sig = kron(P_k_est_sig, P_est_sig);
    end
end

LL_sig = LL_calc(A_obs, P_k_est_sig, psi);


for i=1:k_opt
    if i == 1
        P_k_est = P_est;
    else
        P_k_est = kron(P_k_est, P_est);
    end
end

LL_ori = LL_calc(A_obs, P_k_est, psi);


num_gradient = (LL_sig - LL_ori)./(delta)
approx_gradiet = gradient_calc_v2(P_est, P_k_est, m_opt, k_opt,  G_sample, psi)