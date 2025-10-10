%% Definition of system parameters
B = 5; T = 3; Nt = 16;
data_number  = 1000;
syms x K; semantic_rate = (-0.009854*K^3+0.077*K^2+0.03711*K)/(1+exp(-0.2*x)*exp(-0.1133*K^2+0.6573*K-1.722))+0.01218*K^2-0.1009*K+0.283; 
numerical_semantic_rate = matlabFunction(semantic_rate);
SNR = 0;
sigma = 1/(10^(SNR/10)); 

%% downsample ratio and performance metric
K = 4;
L = 1; J = 1/4^(5-K); M = 1;
a = 0.01218*K^2-0.1009*K+0.283;
b = (-0.009854*K^3+0.077*K^2+0.03711*K)/(exp(-0.1133*K^2+0.6573*K-1.722));
c = 1/(exp(-0.1133*K^2+0.6573*K-1.722));
alpha = 2/log(10);

%% generate channel data
[channel_dataset,precoding_matrix_dataset] = generate_dataset(data_number,Nt,B,T,10,L,J,M,a,b,c,alpha,ones(B,1),sigma);
save(sprintf('./dataset/channel_dataset_Nt_%d_B_%d_T_%d_K_%d_num_%d.mat',Nt,B,T,K,data_number),'channel_dataset','precoding_matrix_dataset');
%load('channel_dataset.mat');