%This is test utility to validate the attack model
clear;
%For reproducibility
%rng(0);
load('syn_1024_1.mat');
alpha = 0;
for alpha = 50:10:100
[G_new, sequence, index_array, max_sequence] = attack_model(G, alpha, 500);
A = adjacency(G);
[N, m] = size(A);
for i = 1:N
    degree_array(i) = sum(A(i,:));
end
%figure(alpha)
degree_victim = degree_array(sequence);
%histogram(degree_victim)
%xlim([0 60])
% plot(degree_victim)
% hold on
plot(max_sequence(1,:))
hold on
plot(max_sequence(2,:))
hold off
end
