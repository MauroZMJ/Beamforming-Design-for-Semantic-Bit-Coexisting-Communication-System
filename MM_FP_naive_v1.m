function V = MM_FP_naive_v1(channel,B,T,L,J,M,a,b,c,alpha,beta,sigma)
Nt = size(channel,1); K = size(channel,2);
bituser_channel = channel(:,1:B); semuser_channel = channel(:,B+1:K); 
[mu,nu,zeta,yeta] = naive_set(bituser_channel,semuser_channel,L,J,M,B,T,sigma,a,b,c,alpha);
[bituser_precoding_normalize,semuser_precoding_normalize] = precoding_construction(bituser_channel,semuser_channel,sigma,Nt,B,T,mu,nu,zeta,yeta,ones(B,1),ones(T,1),1);
[bituser_precoding_po,semuser_precoding_po] = power_optimization(channel,bituser_precoding_normalize,semuser_precoding_normalize,B,T,L,J,M,a,b,c,alpha,beta,sigma);
V = [bituser_precoding_po semuser_precoding_po];
end
function [mu,nu,zeta,yeta] = naive_set(bituser_channel,semuser_channel,L,J,M,B,T,sigma,a,b,c,alpha)
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
    for bit_idx = 1:B
        semuser_matrix = semuser_matrix + yeta(bit_idx)* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    for sem_idx = 1:T
        semuser_matrix = semuser_matrix + nu(sem_idx) * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    semuser_matrix = semuser_matrix + eye_scalar * sigma * eye(Nt);
    reconstruct_semuser_precoding = zeros(Nt,T);
    for sem_idx = 1:T
        reverse_matrix_instance = semuser_matrix + (zeta(sem_idx)-nu(sem_idx))*semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
        precoding_instance = kexi(sem_idx) * reverse_matrix_instance^(-1)* semuser_channel(:,sem_idx);
        reconstruct_semuser_precoding(:,sem_idx) = precoding_instance;
    end
end
function [bituser_precoding_po,semuser_precoding_po] = power_optimization(channel,bituser_precoding,semuser_precoding,B,T,L,J,M,a,b,c,alpha,beta,sigma)
    V_instance = [bituser_precoding semuser_precoding];
    Nt = size(channel,1); K = size(channel,2);
    V = zeros(Nt,K);
    for i = 1:K
       V(:,i) = V_instance(:,i)/norm(V_instance(:,i),2)*exp(-1i*angle(channel(:,i)'*V_instance(:,i))); 
    end
    V = QOS_PowerControl_v4(channel,V,B,T,L,J,M,a,b,c,alpha,beta,sigma);
    bituser_precoding_po = V(:,1:B); semuser_precoding_po = V(:,B+1:B+T);
end
