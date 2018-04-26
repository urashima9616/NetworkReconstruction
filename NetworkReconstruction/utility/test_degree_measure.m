%Step 1: MFNG generation
%Test MFNG_gen
clear;
clc;
N = 1024;
Burnin = 0;
N_samples = 10;
rng(1);
use_synthesized = 1;
use_correct_order = 1;
ini_method = 'forced degree';
forced_degree_lb = 25;
alpha = 20;
time_steps = 10;
target = 2;
%Option 1: generate new graph
if use_synthesized == 0
    [G, P, k_opt, m_opt] = MFNG_gen_v2(N);
    for i=1:k_opt
        if i == 1
            P_k = P;
        else
            P_k = kron(P_k, P);
        end
    end
else
%Option 2: use a synthesized graph

    load('syn_1024_1.mat');
    load('syn_1024_1_P.mat');
    load('syn_1024_1_Pk.mat');
    k_opt = 10;
    m_opt = 2;
    fprintf('Load the network of size %d\n', m_opt^k_opt);
end

D = degree(G);
[B, I] = sort(P_k(target, :));
I = fliplr(I);
ll = zeros(1, N-1);
degree_array = zeros(1, N);
alpha = 20;
%Calc the likelihood
    norm_factor = 0;
    A_obs = adjacency(G);
%Single attack
for i = 1:N
   if A_obs(target, i) == 1
       A_obs(target, i) = 0 ;
       A_obs(i, target) = 0 ;
   end
end


G_test = graph(A_obs);


for degree = 1:N-1
    counter = 0;
    temp = 0;
    A_test = adjacency(G_test);
    norm_factor = 0;
    for kk = 1:length(I)
        if I(kk) == target
            continue
        end
        if counter < degree
            temp = temp + log(P_k(I(kk), target));
            A_test(target, I(kk)) = 1;
            A_test(I(kk), target) = 1; 
        else
            temp = temp + log(1-P_k(I(kk), target));
        end
        counter = counter + 1;
    end
    for i = 1:N
        degree_array(i) = sum(A_test(i,:));
        norm_factor = norm_factor + degree_array(i)^(alpha);
    end
    victim_likelihood = degree_array(target)^(alpha)/norm_factor;
    
    
    ll(1, degree) = temp + log(victim_likelihood);
end
plot(ll(1,1:64))
figure(2)
plot(D)
% 
% fprintf('Attack the network with alpha = %f for %d steps\n', alpha, time_steps);
% %Step 2: Attack the graph
% [G_new, sequence, index_array] = attack_model(G, alpha, time_steps);