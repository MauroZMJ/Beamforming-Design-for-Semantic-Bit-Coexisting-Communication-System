function precoding_matrix = SINR_balance_beamforming(channel,B,T,L,J,M,a,b,c,alpha,beta,sigma,p_sem)
    Nt = size(channel,1); K = size(channel,2);
    precoding_matrix = channel*(channel'*channel)^(-1); precoding_matrix = precoding_matrix ./ sqrt(sum(sum(abs(precoding_matrix).^2)));
    bituser_channel = channel(:,1:B); semuser_channel = channel(:,B+1:K); %seperate the channel of bituser and semuser
    bituser_precoding = precoding_matrix(:,1:B); semuser_precoding = precoding_matrix(:,B+1:B+T);
    semuser_precoding_retain = zeros(size(semuser_precoding));
    for t = 1:T
        semuser_precoding_retain(:,t) =  semuser_precoding(:,t)/norm( semuser_precoding(:,t),2)*(p_sem);
        semuser_precoding(:,t) =  semuser_precoding(:,t)/norm( semuser_precoding(:,t),2)*0;
    end
    end_condition = 0;
    iteration = 0;
    object_performance_list = [];
    cvx_obj_list = [];
    last_performance = 0;
    while end_condition == 0
            [y, z] = update_sinr_term(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T);
            [m, n] = update_fraction_variable(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,y,z,B,T);
            [object_performance,qos,power] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T);
            object_performance = qos;
            [bituser_precoding,semuser_precoding,solving_result,obj_value] = update_precoding_matrix(bituser_channel,semuser_channel,semuser_precoding,sigma,Nt,B,T,y,z,m,n,beta,L,J,M);
            iteration = iteration + 1
            object_performance_list = [object_performance_list,object_performance];
            if iteration == 10
                la = 1;
                %la = for_test(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,gamma,y,z,x,m,n,sigma,a,b,c,alpha,L,J,M,B,T)
                %         plot(sum(object_performance_list,1)); hold on
                %         plot(sum(cvx_obj_list,1))
            end
            if abs(last_performance-sum(object_performance)) <=1e-4 || solving_result == 0 %0 denote solving failure
                end_condition = 1;
            else
                last_performance = sum(object_performance);
            end
    end
    qos
    precoding_matrix = cat(2,bituser_precoding,semuser_precoding_retain);
end
function [y, z] = update_sinr_term(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T)
    y = zeros(B,1);
    z = zeros(B,1);
    %update y
    for b = 1:B
        this_user_channel = bituser_channel(:,b);
        interference = sigma;
        for i = 1:T
            interference = interference + abs(this_user_channel'*semuser_precoding(:,i)).^2;
        end
        for i = 1:B
            if i~=b
                interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
            end
        end
        y(b) = abs(this_user_channel'*bituser_precoding(:,b)).^2/interference;
    end
    %update z
    for b = 1:B
        this_user_channel = bituser_channel(:,b);
        interference = sigma;
        for i = 1:B
            if i~=b
                interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
            end
        end
        z(b) = abs(this_user_channel'*bituser_precoding(:,b)).^2/interference;
    end
end
function [m, n] = update_fraction_variable(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,y,z,B,T)
    m = zeros(B,1);
    n = zeros(B,1);
    %update m
    for b = 1:B
        this_user_channel = bituser_channel(:,b);
        interference = sigma;
        for i = 1:T
            interference = interference + abs(this_user_channel'*semuser_precoding(:,i)).^2;
        end
        for i = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        m(b) = sqrt(1+y(b)).*abs(real(this_user_channel'*bituser_precoding(:,b)))/interference;
    end
    %update n
    for b = 1:B
        this_user_channel = bituser_channel(:,b);
        interference = sigma;
        for i = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        n(b) = sqrt(1+z(b)).*abs(real(this_user_channel'*bituser_precoding(:,b)))/interference;
    end
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
    u = min(rate_list);
    rate_list
    for bit_idx = 1:B
        this_user_channel = bituser_channel(:,bit_idx);
        reconstruct_lamda(bit_idx) = (m(bit_idx)^2*interference_list(bit_idx) + beta(bit_idx)*u-log2(1+y(bit_idx))+y(bit_idx))/(2*m(bit_idx)^2*(1+y(bit_idx))*real(this_user_channel'*inverse_bituser_matrix*this_user_channel));
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
function [bituser_precoding,semuser_precoding,solving_result,obj_value] = update_precoding_matrix(bituser_channel,semuser_channel,semuser_precoding,sigma,Nt,B,T,y,z,m,n,beta,L,J,M)
    reconstruct_lamda = ones(B,1)*1; reconstruct_lamda_old = zeros(B,1);
    reconstruct_nu = 1; 

    while sum(abs(reconstruct_lamda-reconstruct_lamda_old))>1e-6
        reconstruct_lamda_old = reconstruct_lamda;
        [reconstruct_bituser_precoding,inverse_bituser_matrix] = precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,y,z,m,n,reconstruct_lamda,reconstruct_nu,1);    
        reconstruct_lamda = update_lamda(bituser_channel,reconstruct_bituser_precoding,inverse_bituser_matrix,m,y,sigma,beta,B);
    end
    solving_result = 1;
    bituser_precoding = reconstruct_bituser_precoding;
    [~,obj_list,~] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,1,1,1,1,L,J,M,B,T);
    obj_value = sum(obj_list);
end


