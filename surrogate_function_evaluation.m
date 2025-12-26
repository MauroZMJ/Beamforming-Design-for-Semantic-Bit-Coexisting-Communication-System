%In this file, we evaluate the proposed surrogate function. 
clear all
K = 5; %K=4 for another figure.  % downsample times

x = 0.04:0.01:4; 

%parameters of semantic rate function, obtained via data fitting.  
a_list = [0.1495,0.1101,0.0874,0.0747,0.0793,0.1195]; a = a_list(K);
b_list = [0.4659,0.7734,1.34,2.37,3.844,7.33]; b = b_list(K);
c_list = [1.806,1.675,1.971,2.868,4.219,8.297]; c = c_list(K);
alpha_list = [0.1992,0.2198,0.2289,0.2253,0.244,0.2358]; 
alpha = 10*alpha_list(K)/log(10);


%The proposed surrogate function
if alpha <= 1
    x0 = 0.2;
    y1 = a + b./(c+(alpha*x0+(1-alpha)*x)./(x0^(alpha)*x));
    y1_zero = a+b./(c+x0.^(-alpha));
    x0 = 3.0;
    y2 = a + b./(c+(alpha*x0+(1-alpha)*x)./(x0^(alpha)*x));
    y2_zero = a+b./(c+x0.^(-alpha));
else
    x0 = 0.2;
    y1 = a + b./(c+1./(x0^(alpha-1)*(x0+alpha*(x-x0))));
    y1_zero = a+b./(c+x0.^(-alpha));
    x0 = 3.0;
    y2 = a + b./(c+1./(x0^(alpha-1)*(x0+alpha*(x-x0))));
    y2_zero = a+b./(c+x0.^(-alpha));
end

%GT function
y = a+b./(c+x.^(-alpha));


plot(x,y,'-','linewidth',1.5);
hold on

plot(x,y2,'--','linewidth',1.5);
hold on

plot(x,y1,'--','linewidth',1.5);
hold on

plot(0.2,y1_zero,'k-o','MarkerSize',10,'MarkerFaceColor','#EDB120'); hold on
plot(3,y2_zero,'k-o','MarkerSize',10,'MarkerFaceColor','#D95319'); hold on
axis([0,4,0.2,1])
set(gca,'linewidth',1,'fontsize',16,'fontname','Times');
le = legend('Real Semantic Rate','Proposed Surrogate Function($\Gamma_{t_i}^0$=0.2)','Proposed Surrogate Function($\Gamma_{t_i}^0$=3.5)','Tangent Approximation($\Gamma_{t_i}^0$=0.2)','Tangent Approximation($\Gamma_{t_i}^0$=3.5)','fontname','Times','interpreter','latex'); set(le,'linewidth',1,'fontsize',16,'fontname','Times');

grid on; %title('SCA Evaluation','Fontname', '宋体','FontSize',20);
xlabel('$\Gamma_{t_i}$','Fontname', 'Times New Roman','FontSize',20,'interpreter','latex');
ylabel('$\tilde{\epsilon}_K$','Fontname', 'Times New Roman','FontSize',20,'interpreter','latex');
