%This is a test utility to check the quality of the samples drawn 
%Input: attack sequence, Z_samples, original graph G
temp = degree(G);


N = length(temp);
degree_dist = zeros(length(Z_samples)+1, N);


avg_dist = zeros(1, N);
for i = 1: length(Z_samples)
   G_temp = graph(Z_samples{i});
   temp_dist = degree(G_temp);
   %Extract the degree of nodes removed
   degree_dist(i, :) = temp_dist';
   avg_dist = avg_dist + temp_dist';
   %plot(degree_dist(i,:))
   %hold on
end
%figure(1)
%temp = degree(G);
%temp_dist = temp(sequence)';
%plot(temp_dist, 'r')
%hold on

avg_dist = avg_dist./length(Z_samples);
plot(avg_dist)
hold on
plot(temp)


