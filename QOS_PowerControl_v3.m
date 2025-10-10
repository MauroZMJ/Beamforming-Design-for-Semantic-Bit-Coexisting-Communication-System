function precoding_matrix = QOS_PowerControl(channel,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma)
    Nt = size(channel,1); K = size(channel,2);
    bituser_channel = channel(:,1:B); semuser_channel = channel(:,B+1:K); %seperate the channel of bituser and semuser
    hv_matrix = channel'*precoding_matrix;
    end_condition = 0;
    iteration = 0;
    object_performance_list = [];
    cvx_obj_list = [];
    last_performance = 0;
    [p_bit,p_sem] = QOS_PowerControl_feasible_point(channel,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma);
    [~,qos,~] = performance(bituser_channel,semuser_channel,precoding_matrix(:,1:B)*diag(sqrt(p_bit)),precoding_matrix(:,B+1:K)*diag(sqrt(p_sem)),sigma,a,b,c,alpha,L,J,M,B,T);
    if sum(qos>= beta) == 5
        while end_condition == 0
            [gamma,y, z] = update_sinr_term(hv_matrix,p_bit,p_sem,sigma,B,T);
            [x,m, n] = update_fraction_variable(hv_matrix,p_bit,p_sem,sigma,gamma,y,z,B,T,a,b,c,alpha);
            bituser_precoding = precoding_matrix(:,1:B)*diag(sqrt(p_bit));  semuser_precoding = precoding_matrix(:,B+1:K)*diag(sqrt(p_sem));
            [object_performance,qos,power] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T);
            [p_bit,p_sem,solving_result] = update_precoding_matrix_lp(hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,L,J,M,a,b,c,alpha);
            %[p_bit,p_sem,solving_result,obj_value] = update_precoding_matrix_sum(hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,L,J,M,a,b,c,alpha);
            if abs(mean(object_performance)-last_performance) <= 1e-3 && max(abs(qos-beta)) <= 1e-2
                end_condition = 1;
                % [gamma,y, z] = update_sinr_term(hv_matrix,p_bit,p_sem,sigma,B,T);
                % [x,m, n] = update_fraction_variable(hv_matrix,p_bit,p_sem,sigma,gamma,y,z,B,T,a,b,c,alpha);
                % recovery(channel,precoding_matrix,hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,L,J,M,a,b,c,alpha,p_bit,p_sem)
            else
                last_performance = mean(object_performance);
            end
            iteration = iteration + 1;
        end
    else
        p_bit = p_bit; p_sem = ones(T,1)*0;
    end
    bituser_precoding = precoding_matrix(:,1:B)*diag(sqrt(p_bit));  semuser_precoding = precoding_matrix(:,B+1:K)*diag(sqrt(p_sem));
    precoding_matrix(:,1:B) = bituser_precoding; precoding_matrix(:,B+1:K) = semuser_precoding;
end
function [bituser_precoding,semuser_precoding] = recovery(channel,precoding_matrix,hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,L,J,M,a,b,c,alpha,p_bit,p_sem)
    Nt = size(channel,1); K = size(channel,2);
    bituser_channel = channel(:,1:B); semuser_channel = channel(:,B+1:K);
    b_aug = b;
    bituser_precoding = precoding_matrix(:,1:B)*diag(sqrt(p_bit));  semuser_precoding = precoding_matrix(:,B+1:K)*diag(sqrt(p_sem));
    [object_performance,qos,power] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T);
    [D,E,F,G] = cal_scalar(a,b_aug,c,alpha,gamma,T);
    obj_list = zeros(T,1);
    intf = zeros(T,1);
    sinr = zeros(T,1);
    for t = 1:T
        interference = G(t)*sigma+F(t)*p_sem(t)*pow_abs(hv_matrix(B+t,B+t),2);
        for i = 1:T
            if i ~=t
                interference = interference + G(t)*p_sem(i) * pow_abs(hv_matrix(B+t,B+i),2);
            end
        end
        for i = 1:B
            interference = interference + G(t)*p_bit(i)*pow_abs(hv_matrix(B+t,i),2);
        end
        obj_list(t) =  E(t)*(2*x(t)*sqrt(p_sem(t))*real(hv_matrix(B+t,B+t))- x(t)^2 * interference) +D(t);%- o;
        intf(t) = E(t)*p_sem(t)*pow_abs(hv_matrix(B+t,B+t),2)/interference;
    end
end
function [gamma,y, z] = update_sinr_term(hv_matrix,p_bit,p_sem,sigma,B,T)
    gamma = zeros(T,1);
    y = zeros(B,1);
    z = zeros(B,1);
    %update gamma
    for t = 1:T
        interference = sigma;
        for i = 1:T
            if i~=t
                interference = interference + p_sem(i)*abs(hv_matrix(B+t,B+i)).^2;
            end
        end
        for i = 1:B
            interference = interference + p_bit(i)*abs(hv_matrix(B+t,i)).^2;
        end
        gamma(t) = p_sem(t)*abs(hv_matrix(B+t,B+t)).^2/interference;
    end
    %update y
    for b = 1:B
        interference = sigma;
        for i = 1:T
            interference = interference + p_sem(i)*abs(hv_matrix(b,B+i)).^2;
        end
        for i = 1:B
            if i~=b
                interference = interference + p_bit(i)*abs(hv_matrix(b,i)).^2;
            end
        end
        y(b) = p_bit(b)*abs(hv_matrix(b,b)).^2/interference;
    end
    %update z
    for b = 1:B
        interference = sigma;
        for i = 1:B
            if i~=b
                interference = interference + p_bit(i)*abs(hv_matrix(b,i)).^2;
            end
        end
        z(b) = p_bit(b)*abs(hv_matrix(b,b)).^2/interference;
    end
end
function [x,m,n] = update_fraction_variable(hv_matrix,p_bit,p_sem,sigma,gamma,y,z,B,T,a,b,c,alpha)
    x = zeros(T,1);
    m = zeros(B,1);
    n = zeros(B,1);
    [D,E,F,G] = cal_scalar(a,b,c,alpha,gamma,T);
    %update x
    for t = 1:T
        interference = G(t)*sigma+F(t)*p_sem(t)*pow_abs(hv_matrix(B+t,B+t),2);
        for i = 1:T
            if i ~=t
                interference = interference + G(t)*p_sem(i)*abs(hv_matrix(B+t,B+i)).^2; 
            end
        end
        for i = 1:B
            interference = interference + G(t)*p_bit(i)*abs(hv_matrix(B+t,i)).^2;
        end
        x(t) = sqrt(p_sem(t)).*abs(real(hv_matrix(B+t,B+t)))/interference;
    end
    %update m
    for b = 1:B
        interference = sigma;
        for i = 1:T
            interference = interference + p_sem(i)*abs(hv_matrix(b,B+i)).^2;
        end
        for i = 1:B
            interference = interference + p_bit(i)*abs(hv_matrix(b,i)).^2;
        end
        m(b) = sqrt(p_bit(b)*(1+y(b))).*abs(real(hv_matrix(b,b)))/interference;
    end
    %update n
    for b = 1:B
        interference = sigma;
        for i = 1:B
            interference = interference + p_bit(i)*abs(hv_matrix(b,i)).^2;
        end
        n(b) = sqrt(p_bit(b)*(1+z(b))).*abs(real(hv_matrix(b,b)))/interference;
    end
end
function [p_bit,p_sem,solving_result,obj_value] = update_precoding_matrix(hv_matrix,sigma,Nt,B,T,y,z,x,m,n,beta,L,J,M)
    cvx_begin quiet
    variable p_bit(B,1);
    variable p_sem(T,1);
    variable o(1);
    expression obj_list(T,1);
    expression qos_list(B,1);
    for b = 1:B
        interference_m = sigma; interference_n = sigma;
        for i = 1:T
            interference_m = interference_m + p_sem(i)*pow_abs(hv_matrix(b,B+i),2);
        end
        for i = 1:B
            interference_m = interference_m + p_bit(i)*pow_abs(hv_matrix(b,i),2);
            interference_n = interference_n + p_bit(i)*pow_abs(hv_matrix(b,i),2);
        end
        qos_list(b) = J/M*(log2(1+y(b))-y(b)+2*m(b)*sqrt(p_bit(b)*(1+y(b)))*real(hv_matrix(b,b))-m(b)^2*interference_m);
        qos_list(b) = qos_list(b) + (L-J)/M*(log2(1+z(b))-z(b)+2*n(b)*sqrt(p_bit(b)*(1+z(b)))*real(hv_matrix(b,b))-n(b)^2*interference_n);
        qos_list(b) = qos_list(b) - beta(b);
    end

    %calculate the power
    power = 0;
    for i = 1:T
        power = power + p_sem(i);
    end
    for i = 1:B
        power = power + p_bit(i);
    end
    for t = 1:T
        interference = sigma;
        for i = 1:T
            if i ~=t
                interference = interference + p_sem(i) * pow_abs(hv_matrix(B+t,B+i),2);
            end
        end
        for i = 1:B
            interference = interference + p_bit(i)*pow_abs(hv_matrix(B+t,i),2);
        end
        obj_list(t) = 2*x(t)*sqrt(p_sem(t))*real(hv_matrix(B+t,B+t))- x(t)^2 * interference - o;
    end
    maximize(o)
    subject to
    power <= 1
    qos_list >= 0
    obj_list >= 0
    p_sem >= 0
    p_bit >= 0
    cvx_end
    if strcmp(cvx_status,'Solved')==1
        solving_result = 1;
    else
        solving_result = 0;
    end
    obj_value = o;
end
function [D,E,F,G] = cal_scalar(a,b,c,alpha,gamma_0,T)
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
function [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda,mu)
    b_aug = b;
    [D,E,F,G] = cal_scalar(a,b_aug,c,alpha,gamma,T);
    inverse_scalar_list = zeros(B);
    for this_inverse_bit_user = 1:B
        bituser_matrix = 0;
        eye_scalar = 0;
        % inverse matrix calculation
        for bit_idx = 1:B
            scalar = (m(bit_idx)^2*lamda(bit_idx)*J/M + n(bit_idx)^2*lamda(bit_idx)*(L-J))/M;
            eye_scalar = eye_scalar + scalar;
            bituser_matrix = bituser_matrix + scalar* hv_matrix(bit_idx,this_inverse_bit_user)*hv_matrix(bit_idx,this_inverse_bit_user)'; %bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
        end
        for sem_idx = 1:T
            slope = E(sem_idx)*G(sem_idx);
            scalar = slope * x(sem_idx)^2;
            eye_scalar = eye_scalar + scalar;
            bituser_matrix = bituser_matrix + scalar * hv_matrix(B+sem_idx,this_inverse_bit_user)*hv_matrix(B+sem_idx,this_inverse_bit_user)';%semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
        end
        bituser_matrix = bituser_matrix + eye_scalar*sigma;
        inverse_bituser_matrix = bituser_matrix^(-1);
        inverse_scalar_list(this_inverse_bit_user) = inverse_bituser_matrix;
    end
    % bituser precoding recovery
    reconstruct_bituser_precoding = zeros(1,B);
    for bit_idx = 1:B
        inverse_bituser_matrix = inverse_scalar_list(bit_idx);
        scalar = lamda(bit_idx)*m(bit_idx)*sqrt(1+y(bit_idx))*J/M + lamda(bit_idx)*n(bit_idx)*sqrt(1+z(bit_idx))*(L-J)/M;
        precoding_instance = scalar * inverse_bituser_matrix * hv_matrix(bit_idx,bit_idx);%bituser_channel(:,bit_idx);
        reconstruct_bituser_precoding(:,bit_idx) = precoding_instance;
    end

    % inverse matrix calculation
    inverse_scalar_list = zeros(T);
    for this_inverse_sem_user = 1:T
        semuser_matrix = 0;
        %eye_scalar = 0;
        for bit_idx = 1:B
            scalar = lamda(bit_idx)*J/M*m(bit_idx)^2;
            %eye_scalar = eye_scalar + scalar;
            semuser_matrix = semuser_matrix + scalar* hv_matrix(bit_idx,B+this_inverse_sem_user)*hv_matrix(bit_idx,B+this_inverse_sem_user)';
        end
        for sem_idx = 1:T
            slope = E(sem_idx)*G(sem_idx);
            scalar = slope * x(sem_idx)^2;
            %eye_scalar = eye_scalar + scalar;
            semuser_matrix = semuser_matrix + scalar * hv_matrix(B+sem_idx,B+this_inverse_sem_user)*hv_matrix(B+sem_idx,B+this_inverse_sem_user)';
        end
        semuser_matrix = semuser_matrix + eye_scalar * sigma;
        scalar2 = (E(sem_idx)*(G(sem_idx)-F(sem_idx))) * x(sem_idx)^2;
        semuser_matrix = semuser_matrix - scalar2*hv_matrix(B+this_inverse_sem_user,B+this_inverse_sem_user)*hv_matrix(B+this_inverse_sem_user,B+this_inverse_sem_user)';
        inverse_semuser_matrix = semuser_matrix^(-1);
        inverse_scalar_list(this_inverse_sem_user) = inverse_semuser_matrix;
    end

    % semuser precoding recovery
    reconstruct_semuser_precoding = zeros(1,T);
    for sem_idx = 1:T
        inverse_semuser_matrix = inverse_scalar_list(sem_idx);
        coef = max(c*(1-alpha)*gamma(sem_idx).^(alpha)+1,0.01);
        slope = b*alpha*gamma(sem_idx)^(alpha-1)/coef;
        scalar = E(sem_idx) * x(sem_idx); 
        scalar2 = (E(sem_idx)*(G(sem_idx)-F(sem_idx))) * x(sem_idx)^2;
        result_matrix = inverse_semuser_matrix; %+ scalar2/(1-scalar2*semuser_channel(:,sem_idx)'*inverse_semuser_matrix*semuser_channel(:,sem_idx))*inverse_semuser_matrix*semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)'*inverse_semuser_matrix;
        precoding_instance = scalar * result_matrix* hv_matrix(B+sem_idx,B+sem_idx);
        reconstruct_semuser_precoding(:,sem_idx) = precoding_instance;
    end
end
function lamda = update_lambda(u,v,L,J,M,m,n,hv_matrix,bituser_precoding,semuser_precoding,inverse_bituser_matrix,B,T)
    lamda = zeros(B,1);
    for bit_idx = 1:B
        %this_user_channel = bituser_channel(:,bit_idx);
        interference_m = 0; 
        for i = 1:B
            interference_m = interference_m + abs(hv_matrix(bit_idx,i)*bituser_precoding(i)).^2;%abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        interference_n = interference_m; 
        for i = 1:T
            interference_m = interference_m + abs(hv_matrix(bit_idx,B+i)*semuser_precoding(i)).^2;%abs(this_user_channel'*semuser_precoding(:,i)).^2;
        end
        lamda(bit_idx) = (u(bit_idx)+(J)/M*m(bit_idx)^2*interference_m + (L-J)/M*n(bit_idx)^2*interference_n)/(v(bit_idx)*hv_matrix(bit_idx,bit_idx)*bituser_precoding(bit_idx));%real(this_user_channel'*inverse_bituser_matrix*this_user_channel));    
    end
end
function [bituser_precoding,semuser_precoding,solving_result,obj] = update_precoding_matrix_lp(hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,L,J,M,a,b,c,alpha)
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
        %test
        load('./dataset/lambda.mat'); reconstruct_lamda=lambda;p_bit = p_bit'; p_sem=p_sem';%p_bit=sqrt(p_bit); p_sem=sqrt(p_sem);
        [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,reconstruct_lamda,reconstruct_mu);
        snr_scalar = sum(sum(abs(reconstruct_bituser_precoding).^2))+sum(sum(abs(reconstruct_semuser_precoding).^2));
        diff_factor = (snr_scalar-1)*sigma*(J/M*(m.^2)+(L-J)/M*(n.^2));
        reconstruct_lamda = update_lambda(u+diff_factor,v,L,J,M,m,n,hv_matrix,reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,B,T);
        iter = iter + 1;
        if iter > max_iter 
            break;
        end
    end
    [bituser_precoding,semuser_precoding,~,~]=precoding_recovery(hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,reconstruct_lamda,reconstruct_mu);
    solving_result = 1;
    %[obj_list,~,~] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T);
    %obj = sum(obj_list);
    bituser_precoding = real(bituser_precoding).^2/snr_scalar;
    semuser_precoding = real(semuser_precoding).^2/snr_scalar;
    reconstruct_lamda
end

function [p_bit,p_sem,solving_result,obj_value] = update_precoding_matrix_sum(hv_matrix,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,L,J,M,a,b,c,alpha)
    cvx_begin quiet
    variable p_bit(B,1);
    variable p_sem(T,1);
    variable o(1);
    dual variable lambda
    expression obj_list(T,1);
    expression qos_list(B,1);
    b_aug = b;
    [D,E,F,G] = cal_scalar(a,b_aug,c,alpha,gamma,T);
    for b = 1:B
        interference_m = sigma; interference_n = sigma;
        for i = 1:T
            interference_m = interference_m + p_sem(i)^2*pow_abs(hv_matrix(b,B+i),2);
        end
        for i = 1:B
            interference_m = interference_m + p_bit(i)^2*pow_abs(hv_matrix(b,i),2);
            interference_n = interference_n + p_bit(i)^2*pow_abs(hv_matrix(b,i),2);
        end
        qos_list(b) = J/M*(log2(1+y(b))-y(b)+2*m(b)*p_bit(b)*sqrt((1+y(b)))*real(hv_matrix(b,b))-m(b)^2*interference_m);
        qos_list(b) = qos_list(b) + (L-J)/M*(log2(1+z(b))-z(b)+2*n(b)*p_bit(b)*sqrt((1+z(b)))*real(hv_matrix(b,b))-n(b)^2*interference_n);
        qos_list(b) = qos_list(b) - beta(b);
    end
    %calculate the power
    power = 0;
    for i = 1:T
        power = power + p_sem(i)^2;
    end
    for i = 1:B
        power = power + p_bit(i)^2;
    end
    %obj_list = 0;
    for t = 1:T
        interference = G(t)*sigma+F(t)*p_sem(t)^2*pow_abs(hv_matrix(B+t,B+t),2);
        for i = 1:T
            if i ~=t
                interference = interference + G(t)*p_sem(i)^2 * pow_abs(hv_matrix(B+t,B+i),2);
            end
        end
        for i = 1:B
            interference = interference + G(t)*p_bit(i)^2*pow_abs(hv_matrix(B+t,i),2);
        end
        %obj_list(t) =  min(2*x(t)*sqrt(p_sem(t)*E(t))*real(hv_matrix(B+t,B+t))- x(t)^2 * interference +D(t),a+b_aug/c);%- o;
        obj_list(t) =  E(t)*(2*x(t)*(p_sem(t))*real(hv_matrix(B+t,B+t))- x(t)^2 * interference) +D(t);
    end
    %maximize(o)
    maximize(sum(obj_list))
    subject to
    power <= 1
    lambda: qos_list >= 0
    %obj_list >= 0
    p_sem >= 0
    p_bit >= 0
    cvx_end
    if strcmp(cvx_status,'Solved')==1
        solving_result = 1;
    else
        solving_result = 0;
    end
    obj_value = sum(obj_list)
    %obj_list
end