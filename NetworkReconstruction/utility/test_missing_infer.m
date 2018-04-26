clear;
clc;
%Load the observed network
load("fb_graph.mat");

%Attack the network to get the sequence
%G = graph(adjacency_matrix);
adjacency_matrix = full(adjacency(G));
D = degree(G);
header_name = 'facebook';
%User panel
N = length(adjacency_matrix(1,:));
Burnin = 10001;
N_samples = 140000;
rng(3)
beta = 0;

k_opt = 12;
m_opt = 2;
alpha = 100;
loss_percentages = 0.15:0.01:0.15;
P_est_sequence_name = sprintf('%s_est_P_est_sequence%d_burning_%d_samples_beta_%f.mat', header_name, Burnin, N_samples, beta);
psi_name = sprintf('%s_psi_%d_burning_%d_samples_beta_%f.mat', header_name, Burnin, N_samples, beta);
load(P_est_sequence_name);
load(psi_name);

for kk= 1:length(loss_percentages)
    %Attack the network
    loss_percentage = loss_percentages(kk);
    time_steps = floor(loss_percentage*N);
    [G_new, sequence, index_array, max_sequence, prob_sequence] = attack_model_v3(G, alpha, time_steps);
    sequence_missing = zeros(1, time_steps);
    pt = 1;
    for i = 1:N
        if index_array(i) == 0
            sequence_missing(pt) = i;
            pt = pt + 1;
        end
    end

    
    %Get the statistics for sampling
    D_victim = degree(G_new);
    D_effect = D_victim(D_victim>0);
    norm_factor_victim = sum(D_effect.^(alpha));
    sequence_victim = fliplr(sequence);
    A_obs = full(adjacency(G_new));
    index_array_sample = index_array;
    
    %Get the estimated P_est and P_k_est for sampling
    %P_est = P_est_sequence{100,1};
    load('gold_p_est.mat');
    %psi = psi_samples{1,1};
    load('gold_psi.mat')
    
    for i=1:k_opt
        if i == 1
            P_k_est = P_est;
        else
            P_k_est = kron(P_k_est, P_est);
        end
    end
    
    [A_obs_prps, index_array_sample_prps, attack_sequence_prob_prps, flag] = rejection_sampling_attack_opt_nonsym(A_obs, P_k_est, psi, 1, sequence_victim, alpha, index_array_sample, norm_factor_victim, D_victim);
    G_sample_prps = graph(A_obs_prps);
    D_sample_prps = degree(G_sample_prps);
    D_victim_prps = D_sample_prps(sequence_victim);
    D_victim_true = D(sequence_victim);
       

    P_est_base = [0.2863, 0.5795;0.651886, 0.867362]; 
        
    for i=1:k_opt
        if i == 1
            P_k_est_base = P_est_base;
        else
            P_k_est_base = kron(P_k_est_base, P_est_base);
        end
    end
    
   
    
    [A_obs_base, index_array_sample_base, attack_sequence_prob_base, flag] = rejection_sampling_single_v2_random_nonsym(A_obs, P_k_est_base, psi, 1, sequence_missing, alpha, index_array_sample);
    G_sample_base = graph(A_obs_base);
    D_sample_base = degree(G_sample_base);
    D_victim_base = D_sample_base(sequence_victim);
    
    
    
    plot(D_victim_true,'b')
    hold on
    plot(D_victim_prps, 'r')
    hold on
    plot(D_victim_base,'g')
    
    P_est_base = [0.2863, 0.5795;0.651886, 0.867362]; 
    
end