%This is an implementation of the ECM framework 
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
Burnin = 0;
N_samples = 1;
lambda = 10e6;
epsilon = 10e-1;
rng(2);
use_synthesized = 1;
use_correct_order = 1;
ini_method = 'random';
forced_degree_lb = 25;
alpha = 100;
time_steps = 300;
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



global A_obs;
LL = 0;
A_obs = adjacency(G);
batch = 1

while batch < 100
    for index_i = 1:length(P(1,:))
        for index_j = index_i:length(P(1,:))
            iter = 0
            
            while iter < 50
                %E-step:
                %Generate P_k
                gradients = zeros(size(P_est));
                for i=1:k_opt
                    if i == 1
                        P_k_est = P_est;
                    else
                        P_k_est = kron(P_k_est, P_est);
                    end
                end
                
                %Initialize the parameters
                Z_samples = cell(N_samples, 1);
                psi_samples = cell(N_samples, 1);
                attack_sequence_prob_sample = zeros(N_samples, 1);
                LL_sample = zeros(N_samples, 1);
                psi = 1:N;
                
                %Gibbs Sampling:
                
                for kk = 1:Burnin+N_samples
                    LL_sample(kk-Burnin) = LL_calc(A_obs, P_k_est, psi);% + attack_sequence_prob;
                    cur_gradient = gradient_calc(P_est, P_k_est, m_opt, k_opt,  G, psi);
                    gradients = gradients + cur_gradient;
                end
                
                %M_step
                %Construct the Q-function
                gradients = gradients/N_samples;
                P_est_test = P_est ;
                
                P_est_test(index_i, index_j) = P_est_test(index_i, index_j)+ gradients(index_i, index_j)/lambda;
                if P_est_test(index_i, index_j) > 0.9999
                    P_est_test(index_i, index_j) = 0.9999;
                end
                
                if index_i  ~= index_j
                    P_est_test(index_j, index_i) = P_est_test(index_i, index_j);
                end
                P_est = P_est_test;
                
                
                LL_new = sum(LL_sample)/N_samples
                LL_trace = [LL_trace, LL_new];
                
                
                
                %             if abs(LL_new-LL) < epsilon
                %                 break
                %             end
                if LL_new < LL && iter > 25
                    break
                end
                LL = LL_new;
                iter = iter + 1
            end
        end
    end
    batch  = batch + 1
end


