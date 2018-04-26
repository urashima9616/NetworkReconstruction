%This is a utility to calculate the gradient the LL(theta)
%

function gradient_const = gradient_calc_const_v2(P_est, k)

gradient_const = zeros(size(P_est));

const1 = -0.5*k*((sum(sum(P_est)))^(k-1));
const2 = -0.5*k*((sum(sum(P_est.^(2))))^(k-1));
%const1 = -0.5*2*k*((sum(sum(P_est)))^(k-1));
%const2 = -0.5*2*k*((sum(sum(P_est.^(2))))^(k-1));

trace_power_k = -0.5*k*(trace(P_est)^(k-1));
trace_2_power_k = -0.5*k*((trace(P_est.^(2)))^(k-1));

for i = 1:length(P_est(1,:))
    for j = i:length(P_est(1,:))
        if i == j
            gradient_const(i,j) = const1 + trace_power_k + P_est(i,j)*const2 + P_est(i,j)*trace_2_power_k;
        else
            gradient_const(i,j) = 2*const1 + 2*P_est(i,j)*const2;
            gradient_const(j,i) = 2*const1 + 2*P_est(i,j)*const2; 
        end
        
    end
end




