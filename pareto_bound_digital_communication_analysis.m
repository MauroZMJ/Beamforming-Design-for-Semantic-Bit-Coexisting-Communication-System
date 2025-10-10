clear all
%% Definition of system parameters
B = 5; T = 3; Nt = 16;
data_number  = 50;
syms x K; semantic_rate = (-0.009854*K^3+0.077*K^2+0.03711*K)/(1+exp(-0.2*x)*exp(-0.1133*K^2+0.6573*K-1.722))+0.01218*K^2-0.1009*K+0.283; 
numerical_semantic_rate = matlabFunction(semantic_rate);
SNR = 0;
sigma = 1/(10^(SNR/10)); 

%% downsample ratio and performance metric
K = 4;
L = 1; J = 1/4^(5-K); M = 1;
a_list = [0.1495,0.1101,0.0874,0.0747,0.0793,0.1195]; a = a_list(K);
b_list = [0.4659,0.7734,1.34,2.37,3.844,7.33]; b = b_list(K);
c_list = [1.806,1.675,1.971,2.868,4.219,8.297]; c = c_list(K);
alpha_list = [0.1992,0.2198,0.2289,0.2253,0.244,0.2358]; 
alpha = 10*alpha_list(K)/log(10);

%% generate channel data
try
    load(sprintf('./dataset/channel_dataset_Nt_%d_B_%d_T_%d_K_%d_num_%d_SNR_%d.mat',Nt,B,T,K,data_number,SNR));
catch expect
    [channel_dataset,precoding_matrix_dataset] = generate_dataset(data_number,Nt,B,T,10,L,J,M,a,b,c,alpha,ones(B,1),sigma);
    save(sprintf('./dataset/channel_dataset_Nt_%d_B_%d_T_%d_K_%d_num_%d_SNR_%d.mat',Nt,B,T,K,data_number,SNR),'channel_dataset','precoding_matrix_dataset');
end
load('./figures/Digital_Communication_Metric.mat');
%% find the pareto bound
beta_list = 0:0.1:1.5;
performance_list = zeros(size(beta_list));
for idx = 1:length(beta_list)
    beta_now = beta_list(idx)
    beta = beta_now*ones(B,1);
    instance_performance_list = zeros(data_number,1);
   %% Problem solve
    %SCA_FP
    for i = 1:data_number
        i
        channel_instance = squeeze(channel_dataset(i,:,:));
        % initiation
        precoding_matrix = squeeze(precoding_matrix_dataset(i,:,:));
        %precoding_matrix = SINR_balance_beamforming(channel_instance,B,T,L,J,M,a,b,c,alpha,ones(B,1),sigma);
        % calculate the initiate performance
        [object_performance,qos,power] = performance(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,a,b,c,alpha,L,J,M,B,T);
        % problem solve
        if sum(qos>=beta) ==B
            [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP_Maxmin(channel_instance,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma);
        else
            object_performance_list = a * ones(size(object_performance_list));
        end
        [semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = cal_SINR(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,B,T);
        instance_performance_list(i) = mean(Digital_Communication_Metric(10*log10(semuser_sinr)));
    end
    performance_list(idx) = mean(instance_performance_list);
end
%% Performance Evaluation
plot(beta_list,performance_list,'-s','linewidth',1.5,'Color','blue','MarkerSize',10); hold on;
set(gca,'linewidth',1,'fontsize',16,'fontname','Times');
le = legend('Object Performance','SCA Object Performance'); set(le,'linewidth',1,'fontsize',16,'fontname','宋体');
grid on;
%title('卷积编码效果图','Fontname', '宋体','FontSize',20);
xlabel('QOS','Fontname', 'Times New Roman','FontSize',20);
ylabel('Obj Value','Fontname', 'Times New Roman','FontSize',20);
performance_list(i) = object_performance_list(end);