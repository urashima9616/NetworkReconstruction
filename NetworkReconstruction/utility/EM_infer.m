%This is an implementation of the EM framework 
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



%Step 1: MFNG generation
%Test MFNG_gen
clear;
clc;
N = 1024;
Burnin = 100;
N_samples = 200;
rng(1);
use_synthesized = 1;
use_correct_order = 1;
ini_method = 'forced degree';
forced_degree_lb = 25;
alpha = 100;
time_steps = 10;
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

fprintf('Attack the network with alpha = %f for %d steps\n', alpha, time_steps);
%Step 2: Attack the graph
[G_new, sequence, index_array] = attack_model(G, alpha, time_steps);


%Step 3: Network inference
%Initialize the network model parameter matrix
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
        %Do something here to fit the model to the observed graph
        
    otherwise
        fprintf('Illegal initialization method\n');
end





while 1
    %E-step:
    %Generate P_k
    for i=1:k_opt
        if i == 1
            P_k_est = P_est;
        else
            P_k_est = kron(P_k_est, P_est);
        end
    end
   


    %Initialize the  parameters
    Z_samples = cell(N_samples, 1);
    psi_samples = cell(N_samples, 1);
    attack_sequence_prob_sample = zeros(N_samples, 1);
    psi = 1:N;


    %Gibbs Sampling: 

    for kk = 1:Burnin+N_samples
        %Reset G_obs to G_t
        G_obs = G_new;
        index_array_sample = index_array;
        attack_sequence_prob = 0;
        %Part I:  sample P(Z^(tau)|theta(kk), A, G_new, psi^(tau-1))
        %Linear scanning
        if use_correct_order == 0
            for i=1:N
                if index_array_sample(i) == 1
                    continue
                end
                [G_obs, index_array_sample, victim_likelihood] = rejection_sampling(G_obs, P_k_est, psi, i, alpha, index_array_sample);
                attack_sequence_prob = attack_sequence_prob + log(victim_likelihood);
            end    
        else % sanity check. Always infer using the correct ordering
            i = 1;
            while i<length(sequence)
                
                [G_obs, index_array_sample, victim_likelihood, back] = rejection_sampling(G_obs, P_k_est, psi, sequence(i), alpha, index_array_sample);
                if back == 0
                    attack_sequence_prob = attack_sequence_prob + log(victim_likelihood);
                    i  = i + 1;
                else
                    %do some thing here
                end
      
            end   
        end
        
         %###output: Get a G_{t-1:0}^(tau) here !###
        
        
        %Part II:  sample P(psi^(tau)|Z^(tau), theta, A, G_new)
        if use_correct_order == 0
            %psi_sample = Simulate_sampling(G_obs, P_k, index_array);
        else
            psi_sample = psi;
        end
        %###output: Get a psi^(tau) here !###
        
        
        
        psi = psi_sample;
        
        if kk <= Burnin
            continue
        else
            Z_samples{kk-Burnin} = G_obs;
            psi_samples{kk-Burnin} = psi_sample;
            attack_sequence_prob_sample(kk-Burnin) = attack_sequence_prob;
        end

    end

    %M_step
    %Construct the Q-function

    %[P_est, LL_new] = Q_opt(Z_samples, psi_samples, attack_sequence_prob_sample, P_est, k_opt);
    if abs(LL_new-LL) < epsilon
        break
    end
    LL = LL_new;
end
