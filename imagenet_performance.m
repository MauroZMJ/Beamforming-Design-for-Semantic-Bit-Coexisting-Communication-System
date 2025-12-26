%In this file, we provide the results for data fitting. 
load('data/imagenet/performance_index_0.mat'); performance_0 = squeeze(performance(1:end-1,:,2));
load('data/imagenet/performance_index_1.mat'); performance_1 = squeeze(performance(1:end-1,:,2));
load('data/imagenet/performance_index_2.mat'); performance_2 = squeeze(performance(:,:,2));
performance = [performance_0;performance_1;performance_2];
SNR = -30:0.1:30;
%SNR = -30:0.12:30; 
downsample_ratio = 1:1:6;
a = [0.1495,0.1101,0.0874,0.0747,0.0793,0.1195];
b = [0.4659,0.7734,1.34,2.37,3.844,7.33];
c = [1.806,1.675,1.971,2.868,4.219,8.297];
alpha = [0.1992,0.2198,0.2289,0.2253,0.244,0.2358];
figure(1)
SNR_mesh = repmat(SNR',1,length(downsample_ratio)); downsample_ratio_mesh = repmat(downsample_ratio,length(SNR),1);
SNR_mesh_1 = SNR_mesh(:,1); SNR_mesh_2 = SNR_mesh(:,2); SNR_mesh_3 = SNR_mesh(:,3); SNR_mesh_4 = SNR_mesh(:,4); SNR_mesh_5 = SNR_mesh(:,5); SNR_mesh_6 = SNR_mesh(:,6);
SNR_real= 10.^(SNR./10);
SNR_real_mesh = repmat(SNR_real',1,length(downsample_ratio)); downsample_ratio_mesh = repmat(downsample_ratio,length(SNR),1);
SNR_real_mesh_1 = SNR_real_mesh(:,1); SNR_real_mesh_2 = SNR_real_mesh(:,2); SNR_real_mesh_3 = SNR_real_mesh(:,3); SNR_real_mesh_4 = SNR_real_mesh(:,4); SNR_real_mesh_5 = SNR_real_mesh(:,5); SNR_real_mesh_6 = SNR_real_mesh(:,6);

performance_mesh = fliplr(performance(:,:));
performance_mesh_1 = performance_mesh(:,1); performance_mesh_2 = performance_mesh(:,2); performance_mesh_3 = performance_mesh(:,3); performance_mesh_4 = performance_mesh(:,4); performance_mesh_5 = performance_mesh(:,5); performance_mesh_6 = performance_mesh(:,6); 
mesh(SNR_mesh,downsample_ratio_mesh,performance_mesh);
figure(2)
plot(SNR,performance(:,:,1),'linewidth',2);
set(gca,'linewidth',1,'fontsize',16,'fontname','Times');
le = legend('K=1','K=2','K=3','K=4','K=5','K=6'); set(le,'linewidth',1,'fontsize',16,'fontname','Times New Roman');
grid on; %title('SCA Evaluation','Fontname', '宋体','FontSize',20);
xlabel('SNR/dB','Fontname', 'Times New Roman','FontSize',20);
ylabel('SSIM','Fontname', 'Times New Roman','FontSize',20);
%createFit_imagenet(SNR_mesh, downsample_ratio_mesh, performance_mesh);

figure(3)
for i = 1:length(a)
    fitting_result = a(i) + b(i)./(c(i)+SNR_real_mesh_1.^(-10*alpha(i)/log(10)));
    %fitting_result = a(i) + b(i)./(c(i)+exp(-alpha(i).*SNR'));
    plot(SNR',fitting_result-performance(:,7-i,1),'linewidth',1.5); hold on
end
set(gca,'linewidth',1,'fontsize',16,'fontname','Times');
le = legend('K=1','K=2','K=3','K=4','K=5','K=6'); set(le,'linewidth',1,'fontsize',16,'fontname','Times New Roman');
grid on; 
xlabel('SNR','Fontname', 'Times New Roman','FontSize',20);
ylabel('Estimation Error (SSIM)','Fontname', 'Times New Roman','FontSize',20);