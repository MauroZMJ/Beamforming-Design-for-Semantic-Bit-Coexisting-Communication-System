%% Definition of system parameters
B = 5; T = 1; Nt = 16;
data_number  = 10;
beta = 1.2*ones(B,1);
syms x K; semantic_rate = (-0.009854*K^3+0.077*K^2+0.03711*K)/(1+exp(-0.2*x)*exp(-0.1133*K^2+0.6573*K-1.722))+0.01218*K^2-0.1009*K+0.283; 
numerical_semantic_rate = matlabFunction(semantic_rate);
SNR = 0;
sigma = 1/(10^(SNR/10));
L = 1; J = 0.5; M = 1; 
%% generate channel data
channel_dataset = generate_channel(data_number,Nt,B+T,10);
%load('channel_dataset.mat');
%% Problem solve
%WMMSE
%ZF
%MRT 
K = 4;
a = 0.01218*K^2-0.1009*K+0.283;
b = (-0.009854*K^3+0.077*K^2+0.03711*K)/(exp(-0.1133*K^2+0.6573*K-1.722));
c = 1/(exp(-0.1133*K^2+0.6573*K-1.722));
alpha = 2/log(10);
performance_list = zeros(data_number,1);
%SCA_FP
for i = 1:data_number
    channel_instance = squeeze(channel_dataset(i,:,:));
    
    %precoding_matrix =channel_instance*(channel_instance'*channel_instance)^(-1); precoding_matrix = precoding_matrix./norm(precoding_matrix,'fro');
    %precoding_matrix = randn(Nt,B+T) + 1j * randn(Nt,B+T); precoding_matrix = precoding_matrix./norm(precoding_matrix,'fro');
    %load('channel.mat');
    %precoding_matrix =channel_instance*(channel_instance'*channel_instance)^(-1); precoding_matrix = precoding_matrix./norm(precoding_matrix,'fro');
    %% initiation
    %precoding_matrix = ZF_beamforming(channel_instance,B,T,B+T);
    %precoding_matrix = WMMSE_beamforming(channel_instance,B,T,B+T,sigma);
    precoding_matrix = SINR_balance_beamforming(channel_instance,B,T,L,J,M,a,b,c,alpha,beta,sigma);
    
    [object_performance,qos,power] = performance(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,a,b,c,alpha,L,J,M,B,T);
    %[semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = cal_SINR(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,B,T)
    [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP_Maxmin(channel_instance,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma);
    if solving_result == 1
        plot((object_performance_list(:,end)),'-s','linewidth',1.5,'Color','blue','MarkerSize',10); hold on;
        plot((cvx_obj_list(:,end)),'-^','linewidth',1.5,'Color','red','MarkerSize',10);
        set(gca,'linewidth',1,'fontsize',16,'fontname','Times');
        le = legend('Object Performance','SCA Object Performance'); set(le,'linewidth',1,'fontsize',16,'fontname','宋体');
        grid on; 
        %title('卷积编码效果图','Fontname', '宋体','FontSize',20);
        xlabel('Iteration','Fontname', 'Times New Roman','FontSize',20);
        ylabel('Obj Value','Fontname', 'Times New Roman','FontSize',20);
        performance_list(i) = object_performance_list(end);
    end
end
%% Performance Evaluation
performance_list