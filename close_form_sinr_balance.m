clear all
load('workspace.mat')
b = 2.37; %mu=0.6544;
%sigma = 1;
[reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,lamda,nu,0);
reconstruct_lamda = ones(B,1)*0.2; reconstruct_lamda_old = zeros(B,1);
reconstruct_nu = 1; reconstruct_nu_old = 0; 

while sum(abs(reconstruct_lamda-reconstruct_lamda_old))>1e-6
    reconstruct_lamda_old = reconstruct_lamda;
    [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,nu,1);    
    reconstruct_lamda = update_lamda(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B)
end

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
        interference = sigma;
        for b = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,b)).^2;
        end
        interference_list(bit_idx) = interference;
        rate_list(bit_idx) = log2(1+y(bit_idx)) - y(bit_idx) + 2*m(bit_idx)*sqrt(1+y(bit_idx))*real(this_user_channel'*bituser_precoding(:,bit_idx))-m(bit_idx)^2*interference;
    end
    u = max(rate_list);
    rate_list;
    for bit_idx = 1:B
        this_user_channel = bituser_channel(:,bit_idx);
        reconstruct_lamda(bit_idx) = max((m(bit_idx)^2*interference_list(bit_idx) + beta(bit_idx)*u-log2(1+y(bit_idx))+y(bit_idx))/(2*m(bit_idx)^2*(1+y(bit_idx))*real(this_user_channel'*inverse_bituser_matrix*this_user_channel)),0);
    end
end
function [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,lamda,nu,pn)
    bituser_matrix = eye(Nt)*nu;
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