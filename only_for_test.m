% function [precoding_matrix,object_performance_list,cvx_obj_list,solving_result,semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = only_for_test(channel_instance,precoding_matrix_init,B,T,L,K,M,a_list,b_list,c_list,alpha_list,qos,beta_diff,sigma)
% a = a_list(K); b = b_list(K); c = c_list(K); alpha = 10*alpha_list(K)/log(10); J = 1/4^(5-K); 
% %% SCA FP
% if qos(1) >= beta_diff
%     [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP_v2(channel_instance,precoding_matrix_init,B,T,L,J,M,a,b,c,alpha,qos-beta_diff,sigma);
% else
%     object_performance_list = a * ones(T,1);
% end
% [semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = cal_SINR(channel_instance(:,1:B),channel_instance(:,B+1:B+T),precoding_matrix(:,1:B),precoding_matrix(:,B+1:B+T),sigma,B,T);
% end
function [a_list] = only_for_test(channel,precoding_matrix_init,B,T,L,J,M,a,b,c,alpha,beta,sigma,gamma,y,z,x,m,n)
Nt = size(channel,1); K = size(channel,2);
bituser_channel = channel(:,1:B); semuser_channel = channel(:,B+1:K); 
[bituser_precoding,semuser_precoding,solving_result,cvx_obj,lambda] = update_precoding_matrix(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M);
power_scalar = sqrt(sum(sum(abs(bituser_precoding).^2))+sum(sum(abs(semuser_precoding).^2)));
[object_performance,qos,power] = performance(bituser_channel,semuser_channel,bituser_precoding/power_scalar,semuser_precoding/power_scalar,sigma,a,b,c,alpha,L,J,M,B,T);
[mu,nu,zeta,yeta,rho,kexi,eye_scalar] = cal_scalars(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lambda);
[bituser_precoding_new,semuser_precoding_new] = precoding_construction(bituser_channel,semuser_channel,sigma,Nt,B,T,mu,nu,zeta,yeta,rho,kexi,eye_scalar);
tt = sum(mu)+sum(nu);
[bituser_precoding_normalize,semuser_precoding_normalize] = precoding_construction(bituser_channel,semuser_channel,sigma,Nt,B,T,mu/tt,nu/tt,zeta/tt,yeta/tt,ones(B,1),ones(T,1),eye_scalar/tt);
[bituser_precoding_po,semuser_precoding_po] = power_optimization(channel,bituser_precoding,semuser_precoding,B,T,L,J,M,a,b,c,alpha,beta,sigma);
la = 1
[mu,nu,zeta,yeta] = naive_set(bituser_channel,semuser_channel,L,J,M,B,T,x,m,n,sigma,a,b,c,alpha)
[bituser_precoding_normalize,semuser_precoding_normalize] = precoding_construction(bituser_channel,semuser_channel,sigma,Nt,B,T,mu,nu,zeta,yeta,ones(B,1),ones(T,1),1);
[bituser_precoding_po,semuser_precoding_po] = power_optimization(channel,bituser_precoding_normalize,semuser_precoding_normalize,B,T,L,J,M,a,b,c,alpha,beta,sigma);
la =1
%mu_m = mu_m
end
function [bituser_precoding_po,semuser_precoding_po] = power_optimization(channel,bituser_precoding,semuser_precoding,B,T,L,J,M,a,b,c,alpha,beta,sigma)
    V_instance = [bituser_precoding semuser_precoding];
    Nt = size(channel,1); K = size(channel,2);
    V = zeros(Nt,K);
    for i = 1:K
       V(:,i) = V_instance(:,i)/norm(V_instance(:,i),2)*exp(-1i*angle(channel(:,i)'*V_instance(:,i))); 
    end
    V = QOS_PowerControl_v2(channel,V,B,T,L,J,M,a,b,c,alpha,beta,sigma);
    bituser_precoding_po = V(:,1:B); semuser_precoding_po = V(:,B+1:B+T);
end
function [mu,nu,zeta,yeta] = naive_set(bituser_channel,semuser_channel,L,J,M,B,T,x,m,n,sigma,a,b,c,alpha)
    b_aug = b;
    m_eav = zeros(B,1); n_eav = zeros(B,1); x_eav = zeros(T,1); gamma_eav = zeros(T,1);
    %beamforming method
    channel = [bituser_channel semuser_channel];
    %MRT
    V_instance = channel;
    %ZF
    %V_instance = channel*(channel'*channel)^(-1);
    K = size(channel,2); Nt = size(channel,1);
    %normalize
    V = V_instance;
    % for i = 1:K
    %    V(:,i) = V_instance(:,i)/norm(V_instance(:,i),2)*exp(-1i*angle(channel(:,i)'*V_instance(:,i))); 
    % end
    bituser_precoding = V(:,1:B); semuser_precoding = V(:,B+1:B+T);
    for bit_idx = 1:B
        this_user_channel = bituser_channel(:,bit_idx);
        interference_n = sigma; %interference_n = 0;
        for b = 1:B
            interference_n = interference_n + abs(this_user_channel'*bituser_precoding(:,b)).^2;
        end
        interference_m = interference_n;
        for t = 1:T
            interference_m = interference_m + abs(this_user_channel'*semuser_precoding(:,t)).^2;
        end
        m_eav(bit_idx) = real(this_user_channel'*bituser_precoding(:,bit_idx))/interference_m;
        n_eav(bit_idx) = real(this_user_channel'*bituser_precoding(:,bit_idx))/interference_n;
    end

    for sem_idx = 1:T
        this_user_channel = semuser_channel(:,sem_idx);
        interference_x = sigma;
        for b = 1:B
            interference_x = interference_x + abs(this_user_channel'*bituser_precoding(:,b)).^2;
        end
        for t = 1:T
            if t ~= sem_idx
                interference_x = interference_x + abs(this_user_channel'*semuser_precoding(:,t)).^2;
            end
        end
        x_eav(sem_idx) = real(this_user_channel'*semuser_precoding(:,sem_idx))/interference_x;
        gamma_eav(sem_idx) = abs(this_user_channel'*semuser_precoding(:,sem_idx)).^2/interference_x;
    end
    [D,E,F,G ] = cal_weight(a,b_aug,c,alpha,gamma_eav,T);
    for sem_idx = 1:T
        this_user_channel = semuser_channel(:,sem_idx);
        interference_x = G(t)*sigma+F(t)*abs(this_user_channel'*semuser_precoding(:,sem_idx)).^2;
        for b = 1:B
            interference_x = interference_x + G(t)*abs(this_user_channel'*bituser_precoding(:,b)).^2;
        end
        for t = 1:T
            if t ~= sem_idx
                interference_x = interference_x + G(t)*abs(this_user_channel'*semuser_precoding(:,t)).^2;
            end
        end
        x_eav(sem_idx) = real(this_user_channel'*semuser_precoding(:,sem_idx))/interference_x;
    end
    lambda = ones(B,1)*0.2;
    mu = m_eav.^2.*lambda*J/M + n_eav.^2.*lambda*(L-J)/M;
    nu = x_eav.^2.*E.*G;
    zeta = x_eav.^2.*E.*F;
    yeta = m_eav.^2.*lambda*J/M;
    tt = sum(mu) + sum(nu);
    mu = mu/tt; nu = nu/tt; zeta = zeta/tt; yeta = yeta/tt;
end
function [D,E,F,G] = cal_weight(a,b,c,alpha,gamma_0,T)
    if alpha <= 1
        D = a*ones(T,1); E = b*ones(T,1); F = c+(1-alpha)*(gamma_0).^(-alpha); G = alpha*(gamma_0).^(1-alpha);
    else
        slope = max(c * (1-alpha)*gamma_0.^alpha + 1,0.01);
        D = a + b*(1-alpha)*gamma_0.^(alpha)./slope; 
        E = a*b*gamma_0.^(alpha-1)./slope;
        F = a*c*gamma_0.^(alpha-1);
        G = slope;
    end
end
function [mu,nu,zeta,yeta,rho,kexi,eye_scalar] = cal_scalars(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda)
    mu = zeros(B,1); nu = zeros(T,1); zeta = zeros(T,1); yeta = zeros(B,1); rho = zeros(B,1); kexi = zeros(T,1);
    eye_scalar = 0;
    % inverse matrix calculation
    for bit_idx = 1:B
        scalar = (m(bit_idx)^2*lamda(bit_idx)*J/M + n(bit_idx)^2*lamda(bit_idx)*(L-J))/M;
        eye_scalar = eye_scalar + scalar;
        mu(bit_idx) = scalar;
        %bituser_matrix = bituser_matrix + scalar* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    for sem_idx = 1:T
        slope = b*alpha*gamma(sem_idx)^(alpha-1);
        scalar = slope * x(sem_idx)^2;
        eye_scalar = eye_scalar + scalar;
        nu(sem_idx) = scalar;
        %bituser_matrix = bituser_matrix + scalar * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    %bituser_matrix = bituser_matrix + eye_scalar*sigma*eye(Nt);
    %inverse_bituser_matrix = bituser_matrix^(-1);
    % bituser precoding recovery
    %reconstruct_bituser_precoding = zeros(Nt,B);
    for bit_idx = 1:B
        scalar = lamda(bit_idx)*m(bit_idx)*sqrt(1+y(bit_idx))*J/M + lamda(bit_idx)*n(bit_idx)*sqrt(1+z(bit_idx))*(L-J)/M;
        rho(bit_idx) = scalar; 
        %precoding_instance = scalar * inverse_bituser_matrix * bituser_channel(:,bit_idx);
        %reconstruct_bituser_precoding(:,bit_idx) = precoding_instance;
    end

    % inverse matrix calculation
    semuser_matrix = zeros(Nt);
    %eye_scalar = 0;
    for bit_idx = 1:B
        scalar = lamda(bit_idx)*J/M*m(bit_idx)^2;
        yeta(bit_idx) = scalar; 
        %eye_scalar = eye_scalar + scalar;
        %semuser_matrix = semuser_matrix + scalar* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    for sem_idx = 1:T
        %scalar = x(sem_idx)^2*
        coef = max(c*(1-alpha)*gamma(sem_idx).^(alpha)+1,0.01);
        scalar = (b*alpha*gamma(sem_idx)^(alpha-1)*(c*alpha*gamma(sem_idx)^(alpha-1)/coef)) * x(sem_idx)^2;
        zeta(sem_idx) = scalar;
        % slope = b*alpha*gamma(sem_idx)^(alpha-1);
        % scalar = slope * x(sem_idx)^2;
        %eye_scalar = eye_scalar + scalar;
        %semuser_matrix = semuser_matrix + scalar * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    % semuser_matrix = semuser_matrix + eye_scalar * sigma * eye(Nt);
    % inverse_semuser_matrix = semuser_matrix^(-1);
    

    % semuser precoding recovery
    %reconstruct_semuser_precoding = zeros(Nt,T);
    for sem_idx = 1:T
        coef = max(c*(1-alpha)*gamma(sem_idx).^(alpha)+1,0.01);
        slope = b*alpha*gamma(sem_idx)^(alpha-1)/coef;
        scalar = slope * x(sem_idx); 
        kexi(sem_idx) = scalar;
    end
end
function [reconstruct_bituser_precoding,reconstruct_semuser_precoding] = precoding_construction(bituser_channel,semuser_channel,sigma,Nt,B,T,mu,nu,zeta,yeta,rho,kexi,eye_scalar)
    bituser_matrix = zeros(Nt);    
    % inverse matrix calculation
    for bit_idx = 1:B
        bituser_matrix = bituser_matrix + mu(bit_idx)* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    for sem_idx = 1:T
        bituser_matrix = bituser_matrix + nu(sem_idx) * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    bituser_matrix = bituser_matrix + eye_scalar*sigma*eye(Nt);
    inverse_bituser_matrix = bituser_matrix^(-1);
    %bituser precoding recovery
    reconstruct_bituser_precoding = zeros(Nt,B);
    for bit_idx = 1:B
        precoding_instance = rho(bit_idx) * inverse_bituser_matrix * bituser_channel(:,bit_idx);
        reconstruct_bituser_precoding(:,bit_idx) = precoding_instance;
    end

    % inverse matrix calculation
    semuser_matrix = zeros(Nt);
    %eye_scalar = 0;
    for bit_idx = 1:B
        semuser_matrix = semuser_matrix + yeta(bit_idx)* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    for sem_idx = 1:T
        % slope = b*alpha*gamma(sem_idx)^(alpha-1);
        % scalar = slope * x(sem_idx)^2;
        %eye_scalar = eye_scalar + scalar;
        semuser_matrix = semuser_matrix + nu(sem_idx) * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    semuser_matrix = semuser_matrix + eye_scalar * sigma * eye(Nt);
    reconstruct_semuser_precoding = zeros(Nt,T);
    for sem_idx = 1:T
        reverse_matrix_instance = semuser_matrix + (zeta(sem_idx)-nu(sem_idx))*semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
        precoding_instance = kexi(sem_idx) * reverse_matrix_instance^(-1)* semuser_channel(:,sem_idx);
        reconstruct_semuser_precoding(:,sem_idx) = precoding_instance;
    end
    % semuser_matrix = semuser_matrix + eye_scalar * sigma * eye(Nt);
    % inverse_semuser_matrix = semuser_matrix^(-1);
    

    % semuser precoding recovery
    %reconstruct_semuser_precoding = zeros(Nt,T);
    % for sem_idx = 1:T
    %     coef = max(c*(1-alpha)*gamma(sem_idx).^(alpha)+1,0.01);
    %     slope = b*alpha*gamma(sem_idx)^(alpha-1)/coef;
    %     scalar = slope * x(sem_idx); 
    %     scalar2 = (b*alpha*gamma(sem_idx)^(alpha-1)*(1 -c*alpha*gamma(sem_idx)^(alpha-1)/coef)) * x(sem_idx)^2;
    %     result_matrix = inverse_semuser_matrix + scalar2/(1-scalar2*semuser_channel(:,sem_idx)'*inverse_semuser_matrix*semuser_channel(:,sem_idx))*inverse_semuser_matrix*semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)'*inverse_semuser_matrix;
    %     precoding_instance = scalar * result_matrix* semuser_channel(:,sem_idx);
    %     reconstruct_semuser_precoding(:,sem_idx) = precoding_instance;
    % end
end
function [bituser_precoding,semuser_precoding,solving_result,obj,reconstruct_lamda] = update_precoding_matrix(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M)
    u = zeros(B,1); v = zeros(B,1); reconstruct_lamda = ones(B,1)*0.01; reconstruct_lamda_old = zeros(B,1);
    for bit_idx = 1:B
        u(bit_idx) = beta(bit_idx) - (J/M)*(log2(1+y(bit_idx))-y(bit_idx)-m(bit_idx)^2*sigma) - ((L-J)/M)*(log2(1+z(bit_idx))-z(bit_idx)-n(bit_idx)^2*sigma);
        v(bit_idx) = 2*(J/M*m(bit_idx)*sqrt(1+y(bit_idx))+(L-J)/M*n(bit_idx)*sqrt(1+z(bit_idx)))^2;
    end
    reconstruct_mu = 1; reconstruct_mu_old = 0; 
    max_iter = 200; iter = 0;
    %fixed algorithm for the search of lamda
    while sum(abs(reconstruct_lamda-reconstruct_lamda_old))>1e-5
        reconstruct_lamda_old = reconstruct_lamda;
        [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,reconstruct_lamda);
        snr_scalar = sum(sum(abs(reconstruct_bituser_precoding).^2))+sum(sum(abs(reconstruct_semuser_precoding).^2));
        diff_factor = (snr_scalar-1)*sigma*(J/M*(m.^2)+(L-J)/M*(n.^2));
        reconstruct_lamda = update_lambda(u+diff_factor,v,L,J,M,m,n,bituser_channel,semuser_channel,reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,B,T);
        iter = iter + 1;
        if iter > max_iter 
            break;
        end
    end
    [bituser_precoding,semuser_precoding,~,~]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,reconstruct_lamda);
    solving_result = 1;
    [obj_list,~,~] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T);
    obj = sum(obj_list);
    reconstruct_lamda
end
function [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda)
    bituser_matrix = zeros(Nt);
    eye_scalar = 0;
    % inverse matrix calculation
    for bit_idx = 1:B
        scalar = (m(bit_idx)^2*lamda(bit_idx)*J/M + n(bit_idx)^2*lamda(bit_idx)*(L-J))/M;
        eye_scalar = eye_scalar + scalar;
        bituser_matrix = bituser_matrix + scalar* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    for sem_idx = 1:T
        slope = b*alpha*gamma(sem_idx)^(alpha-1);
        scalar = slope * x(sem_idx)^2;
        eye_scalar = eye_scalar + scalar;
        bituser_matrix = bituser_matrix + scalar * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    bituser_matrix = bituser_matrix + eye_scalar*sigma*eye(Nt);
    inverse_bituser_matrix = bituser_matrix^(-1);
    % bituser precoding recovery
    reconstruct_bituser_precoding = zeros(Nt,B);
    for bit_idx = 1:B
        scalar = lamda(bit_idx)*m(bit_idx)*sqrt(1+y(bit_idx))*J/M + lamda(bit_idx)*n(bit_idx)*sqrt(1+z(bit_idx))*(L-J)/M;
        precoding_instance = scalar * inverse_bituser_matrix * bituser_channel(:,bit_idx);
        reconstruct_bituser_precoding(:,bit_idx) = precoding_instance;
    end

    % inverse matrix calculation
    semuser_matrix = zeros(Nt);
    %eye_scalar = 0;
    for bit_idx = 1:B
        scalar = lamda(bit_idx)*J/M*m(bit_idx)^2;
        %eye_scalar = eye_scalar + scalar;
        semuser_matrix = semuser_matrix + scalar* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    for sem_idx = 1:T
        slope = b*alpha*gamma(sem_idx)^(alpha-1);
        scalar = slope * x(sem_idx)^2;
        %eye_scalar = eye_scalar + scalar;
        semuser_matrix = semuser_matrix + scalar * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    semuser_matrix = semuser_matrix + eye_scalar * sigma * eye(Nt);
    inverse_semuser_matrix = semuser_matrix^(-1);

    % semuser precoding recovery
    reconstruct_semuser_precoding = zeros(Nt,T);
    for sem_idx = 1:T
        coef = max(c*(1-alpha)*gamma(sem_idx).^(alpha)+1,0.01);
        slope = b*alpha*gamma(sem_idx)^(alpha-1)/coef;
        scalar = slope * x(sem_idx); 
        scalar2 = (b*alpha*gamma(sem_idx)^(alpha-1)*(1 -c*alpha*gamma(sem_idx)^(alpha-1)/coef)) * x(sem_idx)^2;
        result_matrix = inverse_semuser_matrix + scalar2/(1-scalar2*semuser_channel(:,sem_idx)'*inverse_semuser_matrix*semuser_channel(:,sem_idx))*inverse_semuser_matrix*semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)'*inverse_semuser_matrix;
        precoding_instance = scalar * result_matrix* semuser_channel(:,sem_idx);
        reconstruct_semuser_precoding(:,sem_idx) = precoding_instance;
    end
end
function lamda = update_lambda(u,v,L,J,M,m,n,bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,inverse_bituser_matrix,B,T)
    lamda = zeros(B,1);
    for bit_idx = 1:B
        this_user_channel = bituser_channel(:,bit_idx);
        interference_m = 0; 
        for i = 1:B
            interference_m = interference_m + abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        interference_n = interference_m; 
        for i = 1:T
            interference_m = interference_m + abs(this_user_channel'*semuser_precoding(:,i)).^2;
        end
        lamda(bit_idx) = (u(bit_idx)+(J)/M*m(bit_idx)^2*interference_m + (L-J)/M*n(bit_idx)^2*interference_n)/(v(bit_idx)*real(this_user_channel'*inverse_bituser_matrix*this_user_channel));    
    end
end