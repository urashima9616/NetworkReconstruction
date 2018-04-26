%Draw samples from a pareto distribution
% clear all;
% Xm = 1;
% alpha = 2.1;
% N = 5000;
% Y = randraw('pareto',[Xm, alpha],N);
% Y = sort(Y);
% powers = 1:10:100;
% dist = zeros(length(powers), N);


%Draw samples from a uniform distribution
% clear all;
% Xm = 1;
% alpha = 2.1;
% N = 5000;
% Y = randraw('uniform',[1, 100],N);
% Y = sort(Y);
% powers = 1:10:100;
% dist = zeros(length(powers), N);


%Draw samples from a gaussian distribution
% clear all;
% N = 5000;
% 
% Y = randraw('norm', [25, 4], N);
% Y = sort(Y);
% powers = 1:10:100;
% dist = zeros(length(powers), N);

%Load the samples from the graph
clear;
rng(0);
load('syn_1024_1.mat');
N = 1024;
Y = degree(G);
Y = sort(Y);

powers = 1:10:100;
dist = zeros(length(powers), N);


%Calculate the probability to be target
for i= 1:length(powers)
    power = powers(i);
    Y_aug = Y.^(power);
    norm_factor = sum(Y.^(power));
    dist(i, :) = Y_aug./norm_factor;
    %plot(dist(i,:))
    %pause();
end
