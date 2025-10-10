clear all
system_case = 1;
snr = 10;
if system_case == 1
    dataset_root = './dataset/k_dataset_B_5_T_3_Nt_16/';
else
    dataset_root = './dataset/k_dataset_B_14_T_1_Nt_16/';
end
beta_list = [0.01,0.03,0.05,0.07,0.1];
k_strategy_list = {'random','max','dfo','optimal'};
k_strategy_performance = zeros(length(beta_list),length(k_strategy_list));
for beta_idx = 1:length(beta_list)
    file_name = sprintf('K_dataset_beta_%.3f_snr_%d_.mat',beta_list(beta_idx),snr);
    load(strcat(dataset_root,file_name));
    data_number = length(channel_dataset);
    beta_performance = zeros(data_number,length(k_strategy_list));
    for i = 1:data_number
        k_obj = squeeze(k_obj_dataset(i,1,:));
        k_qos = squeeze(k_qos_dataset(i,1,:));
        for j = 1:length(k_strategy_list)
            switch k_strategy_list{j}
                case 'random'
                    k_predict = randi([1,6]);
                case 'max'
                    k_predict = 6;
                case 'dfo'
                    k_predict = dfo_optimization(4,5,6,k_obj(4),k_obj(5),k_obj(6));
                case 'optimal'
                    [~, k_predict] = max(k_obj);
            end
            beta_performance(i,j) = k_obj(k_predict);
        end
    end
    k_strategy_performance(beta_idx,:) = mean(beta_performance,1);
end
%% Performance Evaluation
plot(beta_list,k_strategy_performance(:,1),'-s','linewidth',1.5,'Color','#0072BD','MarkerSize',10); hold on;
plot(beta_list,k_strategy_performance(:,2),'-d','linewidth',1.5,'Color','#D95319','MarkerSize',10); hold on;
plot(beta_list,k_strategy_performance(:,3),'-^','linewidth',1.5,'Color','#EDB120','MarkerSize',10); hold on;
% plot(SNR_list,SCA_FP_Maxmin_performance_list,'-.','linewidth',1.5,'Color','#7E2F8E','MarkerSize',10); hold on;
% plot(SNR_list,SCA_FP_performance_list,'-*','linewidth',1.5,'Color','#77AC30','MarkerSize',10); hold on;
plot(beta_list,k_strategy_performance(:,4),'-+','linewidth',1.5,'Color','#FFE4E1','MarkerSize',10); hold on;
set(gca,'linewidth',1,'fontsize',16,'fontname','Times');
le = legend('random','max','dfo','optimal'); set(le,'linewidth',1,'fontsize',16,'fontname','宋体');
grid on;
%title('卷积编码效果图','Fontname', '宋体','FontSize',20);
xlabel('SNR','Fontname', 'Times New Roman','FontSize',20);
ylabel('Semantic rate','Fontname', 'Times New Roman','FontSize',20);
function k_predict = dfo_optimization(x1,x2,x3,y1,y2,y3)
    a = ((x1-x2)*(y2-y3)-(y1-y2)*(x2-x3))/((x1-x2)*(x2^2-x3^2)-(x1^2-x2^2)*(x2-x3)); 
    b = ((y1-y2)*(x2^2-x3^2)-(x1^2-x2^2)*(y2-y3))/((x1-x2)*(x2^2-x3^2)-(x1^2-x2^2)*(x2-x3)); 
    c = y1 - a*x1^2 - b*x1;
    if a >=0
        if y1>y3
            k_predict = x1;
        else
            k_predict = x3;
        end
    else
        k_predict = min(max(round(-b/(2*a)),1),6);
    end
end