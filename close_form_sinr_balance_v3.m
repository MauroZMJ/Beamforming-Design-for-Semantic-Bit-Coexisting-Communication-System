clear all
load('workspace.mat')
b = 2.37; %mu=0.6544;
%sigma = 1;
nu = sigma;
[reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,lamda,sigma,1);
reconstruct_lamda = ones(B,1)./B; 
reconstruct_lamda_old = zeros(B,1);
reconstruct_nu = 1; reconstruct_nu_old = 0; 
rate_list = [0;0;0;0;1];

%reconstruct_lamda = [0;0;0;0;1];
% while abs(max(rate_list)-min(rate_list)) > 1e-4%sum(abs(reconstruct_lamda-reconstruct_lamda_old))>1e-6
%     reconstruct_lamda_old = reconstruct_lamda; 
%     %debug
%     %[reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,lamda,nu,0);    
%     %real
%     [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,nu,0);    
%     
%     [rate_list,reconstruct_lamda] = update_lamda_bs(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B,reconstruct_lamda_old);
% end

% while abs(max(rate_list)-min(rate_list)) > 1e-4%sum(abs(reconstruct_lamda-reconstruct_lamda_old))>1e-6
%     reconstruct_lamda_old = reconstruct_lamda; 
%     %debug
%     %[reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,lamda,nu,0);    
%     %real
%     [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,nu,0);    
%     reconstruct_lamda = update_lamda(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B);
%     [rate_list,reconstruct_lamda] = update_lamda_bs(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B,reconstruct_lamda);
% end

while sum(abs(reconstruct_lamda-reconstruct_lamda_old))>1e-6
    reconstruct_lamda_old = reconstruct_lamda; 
    %debug
    %[reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,lamda,nu,0);    
    %real
    [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,nu,0);    
    
    reconstruct_lamda = update_lamda(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B);
    %[rate_list,reconstruct_lamda] = update_lamda_bs(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B,reconstruct_lamda);
end

[reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,nu,0); 
la = 1;
% while sum(abs(reconstruct_nu-reconstruct_nu_old)) >1e-6
%     while sum(abs(reconstruct_lamda-reconstruct_lamda_old))>1e-6
%         reconstruct_lamda_old = reconstruct_lamda;
%         [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,reconstruct_nu,1);    
%         reconstruct_lamda = update_lamda(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B);
%     end
%     reconstruct_lamda_old = zeros(B,1);
%     reconstruct_nu_old = reconstruct_nu;
%     reconstruct_nu = update_nu(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda);
% end
function reconstruct_nu = update_nu(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda)
    nu_low = 0; nu_high = 1;
    [reconstruct_bituser_precoding,~] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,nu_high,0);
    while sum(sum(abs(reconstruct_bituser_precoding).^2))>1
        nu_high = nu_high * 2;
        [reconstruct_bituser_precoding,~] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,nu_high,0);
    end
    while abs(nu_high-nu_low)>1e-4
        nu_avg = (nu_low + nu_high)/2;
        [reconstruct_bituser_precoding,~] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,nu_avg,0);
        if sum(sum(abs(reconstruct_bituser_precoding).^2)) < 1
            nu_high = nu_avg;
        else
            nu_low = nu_avg;
        end
    end
    reconstruct_nu = nu_avg;
end
function reconstruct_lamda = update_lamda(bituser_channel,bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B)
    reconstruct_lamda = zeros(B,1);
    %first calculate u
    rate_list = zeros(B,1);
    interference_list = zeros(B,1);
    for bit_idx = 1:B
        this_user_channel = bituser_channel(:,bit_idx);
        interference = sum(sum(abs(bituser_precoding).^2))*sigma;
        for b = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,b)).^2;
        end
        interference_list(bit_idx) = interference;
        rate_list(bit_idx) = log2(1+y(bit_idx)) - y(bit_idx) + 2*m(bit_idx)*sqrt(1+y(bit_idx))*real(this_user_channel'*bituser_precoding(:,bit_idx))-m(bit_idx)^2*interference;
    end
    %u = 7.7634;
    u = min(rate_list);
    rate_list
    for bit_idx = 1:B
        this_user_channel = bituser_channel(:,bit_idx);
        reconstruct_lamda(bit_idx) = max((m(bit_idx)^2*interference_list(bit_idx) + beta(bit_idx)*u-log2(1+y(bit_idx))+y(bit_idx))/(2*m(bit_idx)^2*(1+y(bit_idx))*real(this_user_channel'*inverse_bituser_matrix*this_user_channel)),0);
    end
    reconstruct_lamda = reconstruct_lamda./sum(reconstruct_lamda);
end
function [rate_list,interference_list] = calculate_rate_qos(bituser_channel,bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B)
    rate_list = zeros(B,1);
    interference_list = zeros(B,1);
    for bit_idx = 1:B
        this_user_channel = bituser_channel(:,bit_idx);
        interference = sum(sum(abs(bituser_precoding).^2))*sigma;
        for b = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,b)).^2;
        end
        interference_list(bit_idx) = interference;
        rate_list(bit_idx) = log2(1+y(bit_idx)) - y(bit_idx) + 2*m(bit_idx)*sqrt(1+y(bit_idx))*real(this_user_channel'*bituser_precoding(:,bit_idx))-m(bit_idx)^2*interference;
    end
end
function [rate_list,reconstruct_lamda] = update_lamda_bs(bituser_channel,bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B,lamda_old)
    reconstruct_lamda = zeros(B,1); 
    for bit_idx = 1:B
        lamda_instance = lamda_old; 
        [rate_list,interference_list] = calculate_rate_qos(bituser_channel,bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B);
        this_user_rate = rate_list(bit_idx);
        if this_user_rate < average(rate_list)
            lamda_user_low = lamda_old(bit_idx); lamda_user_high = 1;
        else
            lamda_user_low = 0;                  lamda_user_high = lamda_old(bit_idx);
        end
        while abs(lamda_user_high-lamda_user_low) > 1e-4
            lamda_user_avg = (lamda_user_low + lamda_user_high)/2;
            lamda_temp = lamda_instance; lamda_temp(bit_idx) = lamda_user_avg;
            [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,bituser_channel(:,1:3),1,16,B,3,y,y,m,m,lamda_temp/sum(lamda_temp),0.01,0);
            [rate_list,interference_list] = calculate_rate_qos(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B);
            this_user_rate = rate_list(bit_idx);
            if this_user_rate < average(rate_list)
                lamda_user_low = lamda_user_avg; 
            else
                lamda_user_high = lamda_user_avg;
            end
        end
        %lamda_old(bit_idx) = lamda_user_avg;
        reconstruct_lamda(bit_idx) = lamda_user_avg;
    end
    rate_list
    reconstruct_lamda = reconstruct_lamda./sum(reconstruct_lamda);
%     %first calculate u
%     [rate_list,interference_list] = calculate_rate_qos(bituser_channel,bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B);
%     %u = 7.7634;
%     u = min(rate_list);
%     rate_list
%     for bit_idx = 1:B
%         this_user_channel = bituser_channel(:,bit_idx);
%         reconstruct_lamda(bit_idx) = max((m(bit_idx)^2*interference_list(bit_idx) + beta(bit_idx)*(sum(rate_list)-rate_list(bit_idx))/(B-1)-log2(1+y(bit_idx))+y(bit_idx))/(2*m(bit_idx)^2*(1+y(bit_idx))*real(this_user_channel'*inverse_bituser_matrix*this_user_channel)),0);
%         %reconstruct_lamda(bit_idx) = max((m(bit_idx)^2*interference_list(bit_idx) + beta(bit_idx)*u-log2(1+y(bit_idx))+y(bit_idx))/(2*m(bit_idx)^2*(1+y(bit_idx))*real(this_user_channel'*inverse_bituser_matrix*this_user_channel)),0);
%     end
%     reconstruct_lamda = reconstruct_lamda./sum(reconstruct_lamda);
end
function [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,lamda,nu,pn)
    bituser_matrix = eye(Nt)*nu*sum(lamda.*m.^2);
    % inverse matrix calculation
    for bit_idx = 1:B
        scalar = m(bit_idx)^2*lamda(bit_idx);
        bituser_matrix = bituser_matrix + scalar* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    inverse_bituser_matrix = bituser_matrix^(-1);
    % bituser precoding recovery
    reconstruct_bituser_precoding = zeros(Nt,B);
    for bit_idx = 1:B
        scalar = lamda(bit_idx)*m(bit_idx)*sqrt(1+y(bit_idx));
        precoding_instance = scalar * inverse_bituser_matrix * bituser_channel(:,bit_idx);
        reconstruct_bituser_precoding(:,bit_idx) = precoding_instance;
    end
    if pn == 1
        reconstruct_bituser_precoding = reconstruct_bituser_precoding./sqrt(sum(sum(abs(reconstruct_bituser_precoding).^2)));
    end
end