 clear;
 rng(1)
 k =10;
 m =2;
 P_est = rand(m);
 P_est = triu(P_est) + triu(P_est)'; %Undirected graph
for i = 1:m
    P_est(i,i)=  P_est(i,i)/2;
end

for i=1:k
    if i == 1
        P_k_est = P_est;
    else
        P_k_est = kron(P_k_est, P_est);
    end
end

%==================Test the likelihood function============================
temp1 = 0;
temp2 = 0;
for i = 1:length(P_k_est(1,:))
    for j = i:length(P_k_est(1,:))
        temp1 = temp1 + log(1-P_k_est(i,j));
        temp2 = -P_k_est(i,j) - 0.5*P_k_est(i,j)^2 + temp2;
    end
end


const1 = 0.5*(   (sum(sum(P_est)))^(k) + trace(P_est)^(k)   );
const2 = 0.5*( (sum(sum((P_est.^(2)))))^(k) +   trace(P_est.^(2))^(k));
temp3 = -(const1 + 0.5*const2);

temp1
temp2
temp3

%==================Test the gradient======================================
%Test the gradient of three methods
%method 1: log(1-P)
gradient_temp1 = zeros(size(P_est));
gradient_temp2 = gradient_temp1;
gradient_temp3 = gradient_temp2;
t_count = zeros(size(P_est));
%Analytical result
for i = 1:length(P_k_est(1,:))
    for j = i:length(P_k_est(1,:))
        t_count = zeros(size(P_est));
        i_labels = calculate_labels(m, k, i);
        j_labels = calculate_labels(m, k, j);
        for kk = 1:length(i_labels)
            t_count(i_labels(kk), j_labels(kk)) =  t_count(i_labels(kk), j_labels(kk)) + 1;
        end
        
        %Update the gradient matrix with the contribution of current edge
        for l = 1:length(gradient_temp1(1,:))
            for n = l : length(gradient_temp1(1,:))
                gradient_temp1(l,n) = gradient_temp1(l,n) + ((P_k_est(i,j)/(1-P_k_est(i,j)))*(-t_count(l,n)))/P_est(l,n);
            end
        end  
    end
end






%Numerical result
P_est_sig = P_est;
delta = 0.001;
P_est_sig(1,2) = P_est_sig(1,2) + delta;
%A_obs = Z_samples{1};

for i=1:k
    if i == 1
        P_k_est_sig = P_est_sig;
    else
        P_k_est_sig = kron(P_k_est_sig, P_est_sig);
    end
end

temp11 = 0;
temp22 = 0;
for i = 1:length(P_k_est_sig(1,:))
    for j = i:length(P_k_est_sig(1,:))
        temp11 = temp11 + log(1-P_k_est_sig(i,j));
        temp22 = -P_k_est_sig(i,j) - 0.5*P_k_est_sig(i,j)^2 + temp22;
    end
end

gradient_temp11 = (temp11-temp1)/delta;
gradient_temp22 = (temp22-temp2)/delta;

const1 = -0.5*k*((sum(sum(P_est)))^(k-1));
const2 = -0.5*k*((sum(sum(P_est.^(2))))^(k-1));
trace_power_k = -0.5*k*(trace(P_est)^(k-1));
trace_2_power_k = -0.5*k*((trace(P_est.^(2)))^(k-1));

for i = 1:length(P_est(1,:))
    for j = i:length(P_est(1,:))
        if i == j
            gradient_temp3(i,j) = const1 + trace_power_k + P_est(i,j)*const2 + P_est(i,j)*trace_2_power_k;
        else
            gradient_temp3(i,j) = const1 + P_est(i,j)*const2;
            gradient_temp3(j,i) = gradient_temp3(i,j); 
        end
        
    end
end

