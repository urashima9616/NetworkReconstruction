%This is a utility to calculate the gradient the LL(theta)
%

function gradient_const = gradient_calc_const_v4(P_est, k)
delta = 0.0001;
gradient_const = zeros(size(P_est));

const1 = 0.5*(   (sum(sum(P_est)))^(k) + trace(P_est)^(k)   );
const2 = 0.5*( (sum(sum((P_est.^(2)))))^(k) +   trace(P_est.^(2))^(k));
LL_orig = -(const1 + 0.5*const2);



for i = 1:length(P_est(1,:))
    for j = i:length(P_est(1,:))
        P_est_sig = P_est;
        P_est_sig(i,j) = P_est_sig(i,j) + delta;
        const1 = 0.5*(   (sum(sum(P_est_sig)))^(k) + trace(P_est_sig)^(k)   );
        const2 = 0.5*( (sum(sum((P_est_sig.^(2)))))^(k) +   trace(P_est_sig.^(2))^(k));
        LL_new = -(const1 + 0.5*const2);
        gradient_const(i,j) = (LL_new-LL_orig)./delta;
        gradient_const(j,i) = gradient_const(i,j);
    end
end




