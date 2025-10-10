clear all
%% Definition of system parameters
B = 14; T = 3; Nt = 16;
data_number  = 50;


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
beta_diff = 0.05;
data_number = 10;
SNR_list = [20]; T_list = [1]; K_list = [1,2,3,4,5];
%SNR_list = [10]; T_list = [1]; K_list = [3,4];
SCA_FP_performance_list = zeros(length(K_list),length(T_list),data_number,length(SNR_list));
sinr_list = zeros(length(K_list),length(T_list),data_number,length(SNR_list));
optimal_k_list = zeros(data_number,length(SNR_list));
for t_idx = 1:length(T_list)
    T = T_list(t_idx)
    channel_dataset = generate_channel(data_number,Nt,B+T,10);
    for idx = 1:length(SNR_list)
        SNR = SNR_list(idx) 
        sigma = 1/(10^(SNR/10));
       %% Problem solve
        for i = 1:data_number
            i
            channel_instance = squeeze(channel_dataset(i,:,:));
            %% SCA_FP initation
            % initiation
            %precoding_matrix = squeeze(precoding_matrix_dataset(i,:,:));
            precoding_matrix_init = SINR_balance_beamforming(channel_instance,B,T,L,0,M,a,b,c,alpha,ones(B,1),sigma,1e-3);
            % calculate the initiate performance
            [object_performance,qos,power] = performance(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix_init(:,1:B),precoding_matrix_init(:,B+1:B+T),sigma,a,b,c,alpha,L,J,M,B,T);
            for K = 1:length(K_list)
                K
                %[precoding_matrix,object_performance_list,cvx_obj_list,solving_result,semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = only_for_test(channel_instance,precoding_matrix_init,B,T,L,K,M,a_list,b_list,c_list,alpha_list,qos,beta_diff,sigma);
                a = a_list(K); b = b_list(K); c = c_list(K); alpha = 10*alpha_list(K)/log(10); J = 1/4^(5-K); 
                %% SCA FP  
                if qos(1) >= beta_diff
                    [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP_v2(channel_instance,precoding_matrix_init,B,T,L,J,M,a,b,c,alpha,qos-beta_diff,sigma);
                else
                    object_performance_list = a * ones(T,1);
                end
                [semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = cal_SINR(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,B,T);
                SCA_FP_performance_list(K,t_idx,i,idx) = mean(object_performance_list(:,end));
                sinr_list(K,t_idx,i,idx) = semuser_sinr;
            end
        end
    end
end
%% Performance Evaluation
optimal_K = zeros(length(T_list),data_number,length(SNR_list));
for t_idx = 1:length(T_list)
    for i = 1:data_number
        for j = 1:length(SNR_list)
            [~,optimal_K(t_idx,i,j)] = max(SCA_FP_performance_list(:,t_idx,i,j));
        end
    end
end
save(sprintf('K_optimization_beta_diff_%.3f.mat',beta_diff),'SCA_FP_performance_list','optimal_K','sinr_list');