clear all
%% Definition of system parameters
B = 5; T = 3; Nt = 16;
data_number  = 1000;
syms x K; semantic_rate = (-0.009854*K^3+0.077*K^2+0.03711*K)/(1+exp(-0.2*x)*exp(-0.1133*K^2+0.6573*K-1.722))+0.01218*K^2-0.1009*K+0.283; 
numerical_semantic_rate = matlabFunction(semantic_rate);


%% downsample ratio and performance metric
K = 4;
L = 1; J = 1/4^(5-K); M = 1;
a_list = [0.1495,0.1101,0.0874,0.0747,0.0793,0.1195]; a = a_list(K);
b_list = [0.4659,0.7734,1.34,2.37,3.844,7.33]; b = b_list(K);
c_list = [1.806,1.675,1.971,2.868,4.219,8.297]; c = c_list(K);
alpha_list = [0.1992,0.2198,0.2289,0.2253,0.244,0.2358]; 
alpha = 10*alpha_list(K)/log(10);

%% generate channel data
channel_dataset = generate_channel(data_number,Nt,B+T,10);
%[channel_dataset,precoding_matrix_dataset] = generate_dataset(data_number,Nt,B,T,10,L,J,M,a,b,c,alpha,ones(B,1),sigma);
%save(sprintf('./dataset/channel_dataset_Nt_%d_B_%d_T_%d_K_%d_num_%d.mat',Nt,B,T,K,data_number),'channel_dataset','precoding_matrix_dataset');
%load(sprintf('./dataset/channel_dataset_Nt_%d_B_%d_T_%d_K_%d_num_%d.mat',Nt,B,T,K,data_number));

%% find the pareto bound
beta_list = 0.5:0.1:1.5;
SNR = 0; 
sigma = 1/(10^(SNR/10));
MRT_performance_list = zeros(size(beta_list));
ZF_performance_list = zeros(size(beta_list));
WMMSE_performance_list = zeros(size(beta_list));
SCA_FP_Maxmin_performance_list = zeros(size(beta_list));
SCA_FP_performance_list = zeros(size(beta_list));
SCA_FP_close_form_performance_list = zeros(size(beta_list));
MM_FP_naive_performance_list = zeros(size(beta_list));
%data_number = 1;
for idx = 1:length(beta_list)
    beta = beta_list(idx)*ones(B,1);
    MRT_instance_performance_list = zeros(data_number,1);
    ZF_instance_performance_list = zeros(data_number,1);
    WMMSE_instance_performance_list = zeros(data_number,1);
    SCA_FP_instance_performance_list = zeros(data_number,1);
    SCA_FP_Maxmin_instance_performance_list = zeros(data_number,1);
    SCA_FP_close_form_instance_performance_list = zeros(data_number,1);
    MM_FP_naive_instance_performance_list = zeros(data_number,1);
   %% Problem solve
    for i = 1:data_number
        i
        channel_instance = squeeze(channel_dataset(i,:,:));
       %% MRT precoding
        precoding_matrix = MRT_beamforming(channel_instance,B,T,L,J,M,a,b,c,alpha,beta,sigma);
        [object_performance,qos,power] = performance(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,a,b,c,alpha,L,J,M,B,T);
        MRT_instance_performance_list(i) = mean(object_performance);    
       %% ZF precoding
        precoding_matrix = ZF_beamforming(channel_instance,B,T,L,J,M,a,b,c,alpha,beta,sigma);
        [object_performance,qos,power] = performance(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,a,b,c,alpha,L,J,M,B,T);
        ZF_instance_performance_list(i) = mean(object_performance);    
       %% WMMSE precoding
        precoding_matrix = WMMSE_beamforming(channel_instance,B,T,L,J,M,a,b,c,alpha,beta,sigma);
        [object_performance,qos,power] = performance(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,a,b,c,alpha,L,J,M,B,T);
        WMMSE_instance_performance_list(i) = mean(object_performance);        
       %% SCA_FP initation
        % initiation
        %precoding_matrix = squeeze(precoding_matrix_dataset(i,:,:));
        precoding_matrix_init = SINR_balance_beamforming_v3(channel_instance,B,T,L,J,M,a,b,c,alpha,ones(B,1),sigma,1e-3);
        % calculate the initiate performance
        [object_performance,qos,power] = performance(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix_init(:,1:B),precoding_matrix_init(:,B+1:B+T),sigma,a,b,c,alpha,L,J,M,B,T);
        %% SCA FP close form 
        tic
        if qos(1) >= beta(1)
            %precoding_matrix = squeeze(precoding_matrix_dataset(i,:,:));
            [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP_v4(channel_instance,precoding_matrix_init,B,T,L,J,M,a,b,c,alpha,beta,sigma);
        else
            object_performance_list = a * ones(size(object_performance));
        end
        toc
        SCA_FP_close_form_instance_performance_list(i) = mean(object_performance_list(:,end));
        
        tic
        precoding_matrix = MM_FP_naive_v1(channel_instance,B,T,L,J,M,a,b,c,alpha,beta,sigma);
        toc
        [object_performance,qos,power] = performance(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,a,b,c,alpha,L,J,M,B,T);  
        MM_FP_naive_instance_performance_list(i) = mean(object_performance); 
        % SCA FP Max Min
%         tic
%         if qos(1) >= beta(1)
%             precoding_matrix = squeeze(precoding_matrix_dataset(i,:,:));
%             [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP_Maxmin(channel_instance,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma);
%         else
%             object_performance_list = a * ones(size(object_performance_list));
%         end
%         toc
%         SCA_FP_Maxmin_instance_performance_list(i) = mean(object_performance_list(:,end));
        
        % SCA FP  
%         tic
%         if qos(1) >= beta(1)
%             precoding_matrix = squeeze(precoding_matrix_dataset(i,:,:));
%             [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP_v2(channel_instance,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma);
%         else
%             object_performance_list = a * ones(size(object_performance_list));
%         end
%         toc
%         SCA_FP_instance_performance_list(i) = mean(object_performance_list(:,end));
    end
    MRT_performance_list(idx) = mean(MRT_instance_performance_list);
    ZF_performance_list(idx) = mean(ZF_instance_performance_list);
    WMMSE_performance_list(idx) = mean(WMMSE_instance_performance_list);
    SCA_FP_Maxmin_performance_list(idx) = mean(SCA_FP_Maxmin_instance_performance_list);
    SCA_FP_performance_list(idx) = mean(SCA_FP_instance_performance_list);
    SCA_FP_close_form_performance_list(idx) = mean(SCA_FP_close_form_instance_performance_list);
    MM_FP_naive_performance_list(idx) = mean(MM_FP_naive_instance_performance_list);
end
%% Performance Evaluation
plot(beta_list,MRT_performance_list,'-s','linewidth',1.5,'Color','#0072BD','MarkerSize',10); hold on;
plot(beta_list,ZF_performance_list,'-d','linewidth',1.5,'Color','#D95319','MarkerSize',10); hold on;
plot(beta_list,WMMSE_performance_list,'-^','linewidth',1.5,'Color','#EDB120','MarkerSize',10); hold on;
plot(beta_list,SCA_FP_Maxmin_performance_list,'-.','linewidth',1.5,'Color','#7E2F8E','MarkerSize',10); hold on;
%plot(beta_list,SCA_FP_performance_list,'-*','linewidth',1.5,'Color','#77AC30','MarkerSize',10); hold on;
plot(beta_list,SCA_FP_close_form_performance_list,'-+','linewidth',1.5,'Color','#FFE4E1','MarkerSize',10); hold on;
plot(beta_list,MM_FP_naive_performance_list,'-*','linewidth',1.5,'Color','#77AC30','MarkerSize',10); hold on;

set(gca,'linewidth',1,'fontsize',16,'fontname','Times');
le = legend('MRT-PC','ZF-PC','WMMSE-PC','SCA-FP-Maxmin','SCA-FP-CF','MM-FP-naive'); set(le,'linewidth',1,'fontsize',16,'fontname','宋体');
grid on;
%title('卷积编码效果图','Fontname', '宋体','FontSize',20);
xlabel('SNR','Fontname', 'Times New Roman','FontSize',20);
ylabel('Semantic rate','Fontname', 'Times New Roman','FontSize',20);