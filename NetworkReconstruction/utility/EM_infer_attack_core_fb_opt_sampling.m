%This is an implementation of the EM framework w/ attack model
%
%Author: Yuankun Xue
%--------------------------------------------------------
%
%
%--Update History:
%12/24/2017:  Originally created by Yuankun Xue

%The EM framework consists of 3 Stages: MFNG generation, MFNG attack, MFNG Inference

%---------------------Step 1: Load the network ----------------------------
clear;
clc;

%load the network:
%G = load('graph.mat');

%Or you can play with graph models

openExample('matlab_featured/BuildWattsStrogatzSmallWorldGraphModelExample')

G = WattsStrogatz(1024,25,0.15);

adjacency_matrix = full(adjacency(G));

%User panel
N = length(adjacency_matrix(1,:));
Burnin = 0;
N_samples = 10000;
N_samples_psi = 10000;
N_batches = 100;
N_opt_steps = 100;
lambda = 10e7;
epsilon = 10e-6;
rng(4);
use_synthesized = 1;
use_correct_order = 1;
use_true_psi = 0;
ini_method = 'random';
forced_degree_lb = 25;
alpha = 100;
beta = 0;
loss_percentage = 0.01;
time_steps = floor(loss_percentage*N);
LL_trace = zeros(1,1);
%Option 1: generate new graph

k_opt = 10;
m_opt = 2;
fprintf('Load the network of size %d\n', m_opt^k_opt);
G = graph(adjacency_matrix);

%Precalculate the labels
labels = zeros(N,k_opt);
for idx = 1:N
    labels(idx, :) = calculate_labels(m_opt, k_opt, idx);
end

%---------------------Step 2: Network Attack------------------------------
fprintf('Attack the network with alpha = %f for %d steps\n', alpha, time_steps);
%Step 2: Attack the graph
[G_new, sequence, index_array] = attack_model_v2(G, alpha, time_steps);
D_victim = degree(G_new);
D_effect = D_victim(D_victim>0);
norm_factor_victim = sum(D_effect.^(alpha));
sequence_victim = fliplr(sequence);


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
        %TODO

    otherwise
        fprintf('Illegal initialization method\n');
end




%>>>>>>>2.EM framework
global A_obs;
global index_array_sample;
LL = 0;
batch = 1;

Z_samples_name =sprintf('graph_samples_attack_%d_loss_%d_nodes_alpha_%d.mat', time_steps, N, alpha);
LL_samples_name =sprintf('LL_samples_attack_%d_loss_%d_nodes_alpha_%d.mat', time_steps, N, alpha);
LL_trace_name = sprintf('LL_trace_attack_%d_loss_%d_nodes_alpha_%d.mat', time_steps, N, alpha);
P_est_sequence_name = sprintf('P_est_sequence_attack_%d_loss_%d_nodes_alpha_%d.mat', time_steps, N, alpha);
P_est_name = sprintf('P_est_attack_%d_loss_%d_nodes_alpha_%d.mat', time_steps, N, alpha);
attack_prob_name = sprintf('attack_prob_attack_%d_loss_%d_nodes_alpha_%d.mat', time_steps, N, alpha);
sequence_name = sprintf('sequence_attack_%d_loss_%d_nodes_alpha_%d.mat', time_steps, N, alpha);
psi_name = sprintf('psi_sequence_attack_%d_loss_%d_nodes_alpha_%d.mat', time_steps, N, alpha);
%Initialize the containers
Z_samples = cell(N_samples, N_batches);
attack_sequence_prob_sample = zeros(N_samples, N_batches);
LL_sample = zeros(N_samples, N_batches);
P_est_sequence = cell(N_opt_steps, N_batches);
psi_samples = cell(N_samples, N_batches);



while batch <= N_batches
    %---------------------------E-step:-----------------------------------
    %Generate P_k

    fprintf('Batch %d\n', batch);
    for i=1:k_opt
        if i == 1
            P_k_est = P_est;
        else
            P_k_est = kron(P_k_est, P_est);
        end
    end


    psi = 1:N;
    sequence_missing = zeros(1, time_steps);
    pt = 1;
    for i = 1:N
        if index_array(i) == 0
            sequence_missing(pt) = i;
            pt = pt + 1;
        end
    end

    %Gibbs Sampling to collect samples:

    for kk = 1:Burnin+N_samples
        %profile on
        %Reset G_obs to G_t
        fprintf('Draw samples.....%d\n', kk);
        A_obs = full(adjacency(G_new));
        index_array_sample = index_array;

       %Part I:  sample P(Z^(tau)|theta(kk), A, G_new, psi^(tau-1))
        if use_correct_order == 0 %Linear scanning
            [A_obs, index_array_sample, attack_sequence_prob, flag] = rejection_sampling_attack_opt(A_obs, P_k_est, psi, 1, sequence_missing, alpha, index_array_sample,  norm_factor_victim, D_victim);

        else % Sanity check. Always infer using the correct ordering
            [A_obs, index_array_sample, attack_sequence_prob, flag] = rejection_sampling_attack_opt(A_obs, P_k_est, psi, 1, sequence_victim, alpha, index_array_sample, norm_factor_victim, D_victim);

        end

        if flag == 0
            fprintf('Failed to generate a sample. Abort.\n');
            return
        end

        %###output: Get a G_{t-1:0}^(tau) here !###


        %Part II:  sample P(psi^(tau)|Z^(tau), theta, A, G_new)
        if use_true_psi == 1
            %psi_sample = Simulate_sampling(G_obs, P_k, index_array);
            psi_sample = psi;
        else
            LL_new = LL_calc(A_obs, P_k_est, psi);
            psi_sample = psi;
            %Take samples on psi
            for s= 1 : N_samples_psi
                %Draw one sample on psi
                [psi_sample, changed, LL_new] = draw_sample_psi(P_k_est, N, LL_new, A_obs, psi_sample);
            end
        end
        %###output: Get a psi^(tau) here !###
        psi = psi_sample;

        if kk <= Burnin
            continue
        else
            Z_samples{kk-Burnin, batch} = sparse(A_obs);
            psi_samples{kk-Burnin,batch} = psi_sample;
            attack_sequence_prob_sample(kk-Burnin, batch) = attack_sequence_prob;
        end
        %profsave
    end

    iter = 0;
    %----------------------------M-step:-----------------------------------
    %Gradient ascent optimization over collected samples
    fprintf('Start to optimize...\n');
    %fprintf('Unoptimized LL value is ')
    %LL_sample(s, batch) = LL_calc(A_obs, P_k_est, psi);
    while iter < N_opt_steps
        %profile on
        gradients = zeros(size(P_est));
        for i=1:k_opt
            if i == 1
                P_k_est = P_est;
            else
                P_k_est = kron(P_k_est, P_est);
            end
        end
        gradient_const = gradient_calc_const_opt(P_est, k_opt, labels);
        for s = 1: N_samples
            A_obs = full(Z_samples{s, batch});
            psi = psi_samples{s,batch};
            LL_sample(s, batch) = LL_calc(A_obs, P_k_est, psi);
            G_sample = graph(A_obs);
            %profile on
            cur_gradient = gradient_calc_opt(P_est, P_k_est, m_opt, k_opt,  G_sample, psi, labels);
            %profsave
            gradients = gradients + cur_gradient;
        end

        %Construct the Q-function
        LL_new = sum(LL_sample(:, batch))/N_samples;
        LL_trace = [LL_trace, LL_new];

        gradients = gradient_const + gradients/N_samples;
        if iter == 0
            momentum = gradients;
            %rms_prop = gradients.^(2);
        else
            momentum = gradients.*(1-beta) + momentum.*(beta);
            %rms_prop = rms_prop.*(beta2) + gradients.^(2)*(1-beta2);
        end

               %P_est_test = P_est + gradients/lambda;
        P_est_test = P_est + momentum/lambda;

        for tt = 1:length(P_est(1,:))
            for ttt = tt:length(P_est(1,:))
                if P_est_test(tt,ttt) > 0.999999
                    P_est_test(tt,ttt) = 0.999999;
                    P_est_test(ttt,tt) = 0.999999;
                end
            end
        end
        P_est = P_est_test;


%         if LL_new < LL && iter >50
%             break
%         end
%            LL = LL_new
            iter = iter + 1  ;

            %profsave
         P_est_sequence{iter, batch} = P_est;
    end
    batch = batch + 1;
end

save(Z_samples_name, 'Z_samples');
save(LL_samples_name, 'LL_sample');
save(LL_trace_name, 'LL_trace');
save(P_est_name, 'P_est');
save(P_est_sequence_name, 'P_est_sequence');
save(attack_prob_name, 'attack_sequence_prob_sample');
save(sequence_name, 'sequence');
save(psi_name, 'psi_samples');
