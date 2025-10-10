clear all
system_case = 1;
snr = 10;
if system_case == 1
    dataset_root = './dataset/k_dataset_B_5_T_3_Nt_16/';
else
    dataset_root = './dataset/k_dataset_B_14_T_1_Nt_16/';
end
a_list = [0.1495,0.1101,0.0874,0.0747,0.0793,0.1195];
b_list = [0.4659,0.7734,1.34,2.37,3.844,7.33];
c_list = [1.806,1.675,1.971,2.868,4.219,8.297];
alpha_list = [0.1992,0.2198,0.2289,0.2253,0.244,0.2358]; 
alpha_list = 10.*alpha_list./log(10);

beta_list = [0.01,0.03,0.05,0.07,0.1];

sinr_list = zeros(length(beta_list),6);
K_list = 1:6;
for beta_idx = 1:length(beta_list)
    file_name = sprintf('K_dataset_beta_%.3f_snr_%d_.mat',beta_list(beta_idx),snr);
    load(strcat(dataset_root,file_name));
    data_number = length(channel_dataset);
    sinr_performance = zeros(data_number,length(beta_list));
    for i = 1:data_number
        k_obj = squeeze(k_obj_dataset(i,1,:));
        k_qos = squeeze(k_qos_dataset(i,1,:));
        for K = 1:length(k_obj)
            sinr_performance(i,K) = (1/(b_list(K)/(k_obj(K)-a_list(K))-c_list(K)))^(1/alpha_list(K));
        end
    end
    sinr_list(beta_idx,:) = mean(sinr_performance,1);
end
%% Performance Evaluation
semilogy(K_list,sinr_list(1,:),'-s','linewidth',1.5,'Color','#0072BD','MarkerSize',10); hold on;
semilogy(K_list,sinr_list(3,:),'-d','linewidth',1.5,'Color','#D95319','MarkerSize',10); hold on;
semilogy(K_list,sinr_list(5,:),'-^','linewidth',1.5,'Color','#EDB120','MarkerSize',10); hold on;

set(gca,'linewidth',1,'fontsize',16,'fontname','Times');
le = legend('random','max','dfo','optimal'); set(le,'linewidth',1,'fontsize',16,'fontname','宋体');
grid on;
%title('卷积编码效果图','Fontname', '宋体','FontSize',20);
xlabel('SNR','Fontname', 'Times New Roman','FontSize',20);
ylabel('Semantic rate','Fontname', 'Times New Roman','FontSize',20);
