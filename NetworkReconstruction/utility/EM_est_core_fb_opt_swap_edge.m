%This is an implementation of the EM framework w/o attack model
%Author: Yuankun Xue
%--------------------------------------------------------
%--Input arguments:
%  num_of_nodes : Number of nodes to be generated.
%  degree_avg: Expected average degree
%  forced_cnnct: The generated network wil be enforced to be connected if
%              non-zero value is set This is achieved by two ways:
%1.Reject disconnected network (forced_cnnct == 1)
%2.Normalize the link probability in a row-wise way as to make sure a node
%  is alway connected to at least one node (forced_cnnct == 2).
%  link_prob: Generating probability measure
%  parts
%
%
%--Update History:
%12/24/2017:  Originally created by Yuankun Xue

%The EM framework consists of 3 Stages: MFNG generation, MFNG attack, MFNG Inference

%---------------------Step 1: Load the network ----------------------------
clear;
clc;

%load the network
load('fb_message.mat')
%load('message_psi_10001_burning_140000_samples_beta_0.000000.mat')
%psi = psi_samples{1,1};
%adjacency_matrix = full(adjacency(G));
%User panel
N = length(adjacency_matrix(1,:));
Burnin = 10001;
N_samples = 100000;
N_batches = 15;
N_opt_steps = 100;
min_step = 0.005;
max_step = 0.1;
lambda = 0.5*10e7;
epsilon = 10e-6;
rng(2);
use_synthesized = 1;
use_correct_order = 1;
ini_method = 'Predefined';
forced_degree_lb = 25;
alpha = 100;
beta = 0;
loss_percentage = 0.02;
time_steps = floor(loss_percentage*N);
LL_trace = zeros(1,1);
%Option 1: generate new graph

P_custom = [0.9, 0.6; 0.6, 0.1];
%P_custom = rand(2);
%P_custom(1,2) = P_custom(2,1);
%P_custom = [0.992310092817853,0.457848428408458;0.457848428408458,0.766268919303346];
%P_custom = [0.984722665853980,0.666646021858268;0.666646021858268,0.408527918448932];
k_opt = 11;
m_opt = 2;
fprintf('Load the network of size %d\n', m_opt^k_opt);
G = graph(adjacency_matrix);
edgelist = table2array(G.Edges); 
len = length(edgelist(:,1));
omega = 0.7;
%Precalculate the labels
labels = zeros(N,k_opt);
for idx = 1:N
    labels(idx, :) = calculate_labels(m_opt, k_opt, idx);
end


%---------------------Step 3: Network Inference----------------------------

%>>>>>>>1.Initialize the network model parameter matrix
fprintf('Network inference starts...\n');
fprintf('Initialize the network model parameter using %s. Expected avg degree is %d\n', ini_method, forced_degree_lb);
switch ini_method
    case 'random'
        %Option1: Random initialization
        P_est = rand(m_opt);
        P_est = triu(P_est) + triu(P_est)'; %Undirected graph
        for i = 1:m_opt
            P_est(i,i)=  P_est(i,i)/2;
        end
    case 'Predefined'
        P_est = P_custom;
    case 'forced degree'
        %Option2: Forced Degree initialization
        ini_avg_degree = 0;
        while ini_avg_degree < forced_degree_lb
            P_est = P + 0.1*rand(m_opt);
            P_est = triu(P_est) + triu(P_est)'; %Undirected graph
            for i = 1:m_opt
                P_est(i,i)=  P_est(i,i)/2;
            end
            for i=1:k_opt
                if i == 1
                    P_k_est = P;
                else
                    P_k_est = kron(P_k_est, P_est);
                end
            end
            ini_avg_degree = sum(sum(P_k_est))/length(P_k_est(1,:));
        end
    case 'optimization'
        %Do something here to fit the model to the observed graph
        
    otherwise
        fprintf('Illegal initialization method\n');
end




%>>>>>>>2.EM framework
LL = 0;
batch = 1;


LL_samples_name =sprintf('message_est_LL_samples_%d_burning_%d_samples_beta_%f.mat', Burnin, N_samples, beta);
LL_trace_name = sprintf('message_est_LL_trace_%d_burning_%d_samples_beta_%f.mat',  Burnin, N_samples, beta);
P_est_sequence_name = sprintf('message_est_P_est_sequence%d_burning_%d_samples_beta_%f.mat', Burnin, N_samples, beta);
P_est_name = sprintf('message_P_est_%d_burning_%d_samples_beta_%f.mat', Burnin, N_samples, beta);
psi_name = sprintf('message_psi_%d_burning_%d_samples_beta_%f.mat', Burnin, N_samples, beta);
%Initialize the containers
LL_sample = zeros(N_samples, N_batches);
P_est_sequence = cell(N_opt_steps, N_batches);
psi_samples = cell(N_opt_steps, N_batches);




psi = 1:N;

iter = 0;
%----------------------------M-step:-----------------------------------
%Gradient ascent optimization over collected samples
fprintf('Start to optimize...\n');
%LL_sample(s, batch) = LL_calc(A_obs, P_k_est, psi);
while iter < N_opt_steps
%     profile on
    tic
    %profile on
    accept_cnt = 0;
    gradients = zeros(size(P_est));
    for i=1:k_opt
        if i == 1
            P_k_est = P_est;
        else
            P_k_est = kron(P_k_est, P_est);
        end
    end
%     if mod(iter, 100) == 0
%         fprintf('Likelihood at %d iterations is....\n', iter);
%         LL_new = LL_calc(adjacency(G), P_k_est, psi)
%     end
    A_obs = full(adjacency(G));
    LL_new = LL_calc(A_obs, P_k_est, psi)
    LL_trace = [LL_trace; LL_new];
    
    gradient_const = gradient_calc_const_opt(P_est, k_opt, labels);
    %gradient_const = gradient_calc_const_opt_nonsym(P_est, k_opt, labels);
    %Take samples on psi
    tic
    for s= 1 : N_samples+ Burnin
        %profile on
        if s >= 1 + Burnin
            [psi_new, changed, LL_new, idx_j, idx_k] = draw_sample_psi_combined(P_k_est, N, LL_new, A_obs, psi, edgelist, len, omega);

            if s == 1 + Burnin
                 cur_gradient = gradient_calc_opt(P_est, P_k_est, m_opt, k_opt,  G, psi_new, labels);
                 %cur_gradient = gradient_calc_opt_nonsym(P_est, P_k_est, m_opt, k_opt,  G, psi_new, labels);
                 %Debug use
                 num_gradient = num_grad_calc(A_obs, P_est, P_k_est, psi_new, m_opt, k_opt);
            else
                if changed
                    accept_cnt = accept_cnt + 1;
                    %Need to incrementally update the gradient
                    cur_gradient = update_gradient(N, A_obs, P_est, P_k_est, labels, cur_gradient, idx_j, idx_k, psi, psi_new);
                    %Debug use
                    test_gradient = gradient_calc_opt(P_est, P_k_est, m_opt, k_opt,  G, psi_new, labels) + gradient_const;
                    num_gradient = num_grad_calc(A_obs, P_est, P_k_est, psi_new, m_opt, k_opt)
                    test_gradient + gradient_const
                    cur_gradient + gradient_const
                    num_gradient
                    
                end
                    
            end
            
            
%             
%             if changed
%                 %Need to incrementally update the gradient
%                 new_cur_grad = update_gradient(N, A_obs, P_est, P_k_est, labels, cur_gradient, idx_j, idx_k, psi, psi_new);
% 
%                 cur_gradient = gradient_calc_opt(P_est, P_k_est, m_opt, k_opt,  G, psi_new, labels);
%             elseif s == 1 + Burnin
%                 cur_gradient = gradient_calc_opt(P_est, P_k_est, m_opt, k_opt,  G, psi, labels);
%             end
%             
            psi = psi_new;
                
        %profsave
            gradients = gradients + cur_gradient;
        else
           %Draw one sample on psi
           [psi, changed, LL_new] = draw_sample_psi_combined(P_k_est, N, LL_new, A_obs, psi, edgelist, len, omega);

           %LL_new
        end
    end
    toc
    psi_samples{+1, 1} = psi;
    accept_rate = accept_cnt./N_samples;
    
    gradients = gradient_const + (gradients./N_samples);
    if iter == 0
        momentum = gradients;
        %rms_prop = gradients.^(2);
    else
        momentum = gradients.*(1-beta) + momentum.*(beta);
        %rms_prop = rms_prop.*(beta2) + gradients.^(2)*(1-beta2);
    end
    
    
    %P_est_test = P_est + gradients/lambda;
    if iter < 100
        lambda = 10e6;
    else
        lambda = 0.5*10e7;
    end
    delta_P = momentum/lambda;
    for g_i = 1: length(P_est(1,:))
        for g_j = 1:length(P_est(1,:))
            if abs(delta_P(g_i, g_j)) < min_step
                if delta_P(g_i, g_j) > 0
                    delta_P(g_i, g_j) = min_step;
                else
                    delta_P(g_i, g_j) = -min_step;
                end
            elseif abs(delta_P(g_i, g_j)) > max_step
                if delta_P(g_i, g_j) > 0
                    delta_P(g_i, g_j) = max_step;
                else
                    delta_P(g_i, g_j) = -max_step;
                end
            end
            
        end
    end
    %P_est_test = P_est + momentum/lambda;
    
    P_est_test = P_est + delta_P;


    
    for tt = 1:length(P_est(1,:))
        for ttt = tt:length(P_est(1,:))
            if P_est_test(tt,ttt) > 0.999999
                P_est_test(tt,ttt) = 0.999999;
                P_est_test(ttt,tt) = 0.999999;
            end
        end
    end
    fprintf("Swap acceptance rate %f\n", accept_rate);
    fprintf("(1,1): %f + %f -> %f\t Gradient: %f \n", P_est(1,1),delta_P(1,1), P_est_test(1,1), momentum(1,1));
    fprintf("(1,2): %f + %f -> %f\t Gradient: %f \n", P_est(1,2),delta_P(1,2), P_est_test(1,2), momentum(1,2));
    fprintf("(2,2): %f + %f -> %f\t Gradient: %f \n", P_est(2,2),delta_P(2,2), P_est_test(2,2), momentum(2,2));
    P_est = P_est_test;
    
    
    
    %            LL = LL_new
    iter = iter + 1
    P_est_sequence{iter, batch} = P_est;
    toc
%     profsave
end

save(P_est_name, 'P_est');
save(P_est_sequence_name, 'P_est_sequence');
save(psi_name, 'psi_samples');



