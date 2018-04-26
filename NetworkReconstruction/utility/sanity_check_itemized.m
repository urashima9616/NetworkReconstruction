%Gradient sanity check on the part dependent on the adjmatrix
%clear;
%load('test_P_est.mat')
load('test_adjmatrix.mat')
G_sample = graph(A_obs);
m = 2;
k = 10;
%rng(6);
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
%=====================Gradient method 2 based on analytical result=========

%Get the edgelist from the graph
e = G_sample.Edges;
edgelist = table2array(e);%
LL_orig = 0;
LL_delta = 0;
psi = 1:1024;
gradient = zeros(size(P_est));
%Iterate on the edge list
for i = 1:length(edgelist(:,1))
    %Get the edge name
    %edge_name = sprintf('$d%d', edgelist(i,1), edgelist(i,2));
    %tf = isKey(mapObj, edge_name);
    
    %If it exists
    mapped_i = psi( edgelist(i,1));
    mapped_j = psi( edgelist(i,2));
    
    LL_orig = LL_orig + log(P_k_est(mapped_i, mapped_j)) - log(1-P_k_est(mapped_i, mapped_j));
    
    
    i_labels = calculate_labels(m, k, mapped_i);
    j_labels = calculate_labels(m, k, mapped_j);
    
%     if tf == 0
%         mapObj(edge_name) = zeros(size(P_est));
%     end
    
    temp_cnt = zeros(size(gradient));
    
    %Get the profile of label hittings 
    for kk = 1:length(i_labels)
        temp_cnt(i_labels(kk), j_labels(kk)) =  temp_cnt(i_labels(kk), j_labels(kk)) + 1;
    end
    %mapObj(edge_name) = temp_cnt;
    
    %Update the gradient matrix with the contribution of current edge
    for l = 1:length(gradient(1,:))
        for n = l : length(gradient(1,:))
            %(P_uv/1-P_uv + 1)*(t/P_i,j)
            gradient(l,n) = gradient(l,n) + (temp_cnt(l,n)/P_est(l,n))*((1-P_k_est(mapped_i, mapped_j))^(-1));
            gradient(n,l) = gradient(l,n);
            %gradient(l,n) = gradient(l,n) - 2*(temp_cnt(l,n)/P_est(l,n));
        end
    end
end

%Numerical result
gradient_num = zeros(size(P_est));
%Calculate the numerical gradient of const term
gradient_const_num = zeros(size(P_est));
LL_const = 0;
delta = 0.001;

for i = 1:length(P_k_est(1,:))
    for j = i:length(P_k_est(1,:))
        LL_const = LL_const + log(1-P_k_est(i,j));
    end
end





for tt = 1:length(P_est)
    for ttt = tt:length(P_est)
        LL_delta = 0;
        LL_const_delta = 0;
        P_est_sig = P_est;
        P_est_sig(tt,ttt) = P_est_sig(tt,ttt) + delta;
        %A_obs = Z_samples{1};
        
        %Generate the link measure
        for i=1:k
            if i == 1
                P_k_est_sig = P_est_sig;
            else
                P_k_est_sig = kron(P_k_est_sig, P_est_sig);
            end
        end
        %Calculate the LL value
        for i = 1:length(P_k_est_sig(1,:))
            for j = i:length(P_k_est_sig(1,:))
                LL_const_delta = LL_const_delta + log(1-P_k_est_sig(i,j));
            end
        end
        
        gradient_const_num(tt,ttt) = (LL_const_delta-LL_const)./delta;
        
        
        
        for i = 1:length(edgelist(:,1))
            mapped_i = psi( edgelist(i,1));
            mapped_j = psi( edgelist(i,2));
            
            LL_delta = LL_delta + log(P_k_est_sig(mapped_i, mapped_j)) - log(1-P_k_est_sig(mapped_i, mapped_j));
        end
        
        %Numerical gradient of the edge term
        gradient_num(tt,ttt) = (LL_delta - LL_orig)/delta;
    end
end



%Verify the rewriting of log-likelihood function
LL_test = LL_const + LL_orig;
LL_all = LL_calc(A_obs, P_k_est, psi);
gradient_const = gradient_calc_const_v4(P_est, k);

gradient_final_approx = gradient + gradient_const;
gradient_final_num = gradient_num + gradient_const_num;




