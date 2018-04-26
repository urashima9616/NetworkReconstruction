%Test MFNG_gen
clear;
clc;
N = 1024;
[G, P, k_opt, m_opt] = MFNG_gen_v2(N);
for i=1:k_opt
    if i == 1
        P_k = P;
    else
        P_k = kron(P_k, P);
    end
end
D = degree(G);
plot(G);
