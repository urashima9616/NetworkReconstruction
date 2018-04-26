%This is a test utility to check the quality of the samples drawn 
%Input: attack sequence, Z_samples, original graph G
N = length(sequence);
degree_dist = zeros(length(Z_samples(:,4))+1, N);
avg_dist = zeros(1, N);


for i = 1: length(Z_samples)
   G_temp = graph(Z_samples{i,4});
   temp = degree(G_temp);
   %Extract the degree of nodes removed
   degree_dist(i, :) = temp(sequence)';
   avg_dist = avg_dist + temp(sequence)';
   %plot(degree_dist(i,:))
   %hold on
end
%figure(1)
avg_dist = avg_dist./length(Z_samples(:,9));
plot(avg_dist)
hold on

 temp = degree(G);
 temp_dist = temp(sequence)';
 plot(temp_dist, 'r')
 hold on




