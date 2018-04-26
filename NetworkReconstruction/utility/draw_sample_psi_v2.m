function [new_psi, changed, LL_new, idx_j, idx_k] = draw_sample_psi_v2(P_k_est, N, LL_prev, A_obs, psi)
%Step 1: draw two sample of index
changed = 0;
idx_j = floor(rand(1)*N + 1);
idx_k = floor(rand(1)*N + 1);
while idx_k == idx_j
    idx_k = floor(rand(1)*N + 1);
end
new_psi = psi;
temp = new_psi(1, idx_j);
new_psi(1, idx_j) = new_psi(1, idx_k);
new_psi(1, idx_k) = temp;

LL_dice = log(rand(1));
%Calculate the new ll
LL_new = LL_prev;
for i=1:N
    if A_obs(i, idx_j) == 1
        diff = - log(P_k_est(psi(i), psi(idx_j) )) + log(P_k_est(new_psi(i), new_psi(idx_j) ));
        LL_new = LL_new + diff;
    elseif A_obs(i, idx_j) == 0 && i ~= idx_j
        diff = - log(1-P_k_est(psi(i), psi(idx_j) )) + log(1-P_k_est(new_psi(i), new_psi(idx_j) ));
        LL_new = LL_new + diff;
    end
    
    if A_obs(i, idx_k) == 1
        diff = - log(P_k_est(psi(i), psi(idx_k) )) + log(P_k_est(new_psi(i), new_psi(idx_k) ));
        LL_new = LL_new + diff;
    elseif A_obs(i, idx_k) == 0 && i ~= idx_k
        diff = - log(1-P_k_est(psi(i), psi(idx_k) )) + log(1-P_k_est(new_psi(i), new_psi(idx_k) ));
        LL_new = LL_new + diff;
    end
end

if A_obs(idx_j, idx_k) == 1
    LL_new = LL_new + log(P_k_est(psi(idx_k), psi(idx_j))) - log(P_k_est(new_psi(idx_k), new_psi(idx_j)));
else
    LL_new = LL_new + log(1-P_k_est(psi(idx_k), psi(idx_j))) - log(1-P_k_est(new_psi(idx_k), new_psi(idx_j)));
end


if LL_dice < LL_new-LL_prev
    changed = 1;
else
    new_psi = psi;
    LL_new = LL_prev;
end