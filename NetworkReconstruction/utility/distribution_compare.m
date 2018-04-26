clear;
clc;
%User panel
N = 1024;
Burnin = 0;
N_samples = 10;
N_batches = 15;
N_opt_steps = 50;
lambda = 0.5*10e6;
epsilon = 10e-6;
rng(1);
use_synthesized = 1;
use_correct_order = 1;
ini_method = 'random';
forced_degree_lb = 25;
alpha = 100;
beta = 0.9;
loss_percentage = 0.2;
time_steps = floor(loss_percentage*N);
LL_trace = zeros(1,1);
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

D1 = degree(G);
[f1, x1] = ksdensity(D1);
plot(x1, f1, 'DisplayName','Original Network')
hold on
alphas = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10];
for i = 1:length(alphas)
    alpha = alphas(i);
header = sprintf("alpha = %d", alpha);
%---------------------Step 2: Network Attack------------------------------
fprintf('Attack the network with alpha = %f for %d steps\n', alpha, time_steps);
%Step 2: Attack the graph
[G_new, sequence, index_array] = attack_model_v2(G, alpha, time_steps);
D2 = D1(sequence);
[f2, x2] = ksdensity(D2);
plot(x2, f2, 'DisplayName', header)
hold on
end
legend('show')






