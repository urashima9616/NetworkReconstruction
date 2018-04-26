%Test MFNG_gen
N = 1024;
avg_d = 4;
G = MFNG_gen(N, avg_d, 0);
D = degree(G);
plot(G);
