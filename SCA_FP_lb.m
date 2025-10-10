function [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP_lb(channel,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma)
Nt = size(channel,1); K = size(channel,2);
bituser_channel = channel(:,1:B); semuser_channel = channel(:,B+1:K); %seperate the channel of bituser and semuser
bituser_precoding = precoding_matrix(:,1:B); semuser_precoding = precoding_matrix(:,B+1:K);
end_condition = 0;
iteration = 0;
object_performance_list = [];
cvx_obj_list = [];
last_performance = 0;
max_iteration = 100;
while end_condition == 0
    [gamma, y, z] = update_sinr_term(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T);
    [x, m, n] = update_fraction_variable(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,gamma,y,z,a,b,c,alpha,B,T);
    power_scalar = sqrt(sum(sum(abs(bituser_precoding).^2))+sum(sum(abs(semuser_precoding).^2)));
    [object_performance,qos,power] = performance(bituser_channel,semuser_channel,bituser_precoding/power_scalar,semuser_precoding/power_scalar,sigma,a,b,c,alpha,L,J,M,B,T);
    %[bit_qos,sem_sinr] = only_for_test(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M);
    %cvx_obj = cal_cvx_obj(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M);
    [bituser_precoding,semuser_precoding,solving_result,cvx_obj] = update_precoding_matrix(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M);
    %[semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = cal_SINR(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T)
    iteration = iteration + 1
    object_performance_list = [object_performance_list,object_performance];
    cvx_obj_list = [cvx_obj_list,cvx_obj];
    if iteration == 4
        la = 1;
        %plot(sum(object_performance_list,1)); hold on
        %plot(sum(cvx_obj_list,1))
    end
    if abs(last_performance-sum(object_performance)) <=1e-4 || solving_result == 0 || iteration>= max_iteration %0 denote solving failure
        end_condition = 1;
    else
        last_performance = sum(object_performance);
    end
end
power_scalar = sqrt(sum(sum(abs(bituser_precoding).^2))+sum(sum(abs(semuser_precoding).^2)));
precoding_matrix(:,1:B) = bituser_precoding/power_scalar;  precoding_matrix(:,B+1:K) = semuser_precoding/power_scalar;
end
function gamma = update_gamma(~,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T)
    gamma = zeros(T,1);
        %update gamma
    for t = 1:T
        this_user_channel = semuser_channel(:,t);
        interference = sigma;
        for i = 1:T
            if i~=t
                interference = interference + abs(this_user_channel'*semuser_precoding(:,i)).^2;
            end
        end
        for i = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        gamma(t) = abs(this_user_channel'*semuser_precoding(:,t)).^2/interference;
    end
end
function x = update_x(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T)
    x = zeros(T,1); 
    %update x
    for t = 1:T
        this_user_channel = semuser_channel(:,t);
        interference = sigma;
        for i = 1:T
            if i~=t
                interference = interference + abs(this_user_channel'*semuser_precoding(:,i)).^2;
            end
        end
        for i = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        x(t) = abs(real(this_user_channel'*semuser_precoding(:,t))) / interference;
    end
end
function [y,z] = update_y_z(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T)
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
function [m,n] = update_m_n(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,y,z,B,T)
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
            if i~=b
                interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
            end
        end
        m(b) = sqrt(1+y(b)).*abs(real(this_user_channel'*bituser_precoding(:,b)))/interference;
    end
    %update n
    for b = 1:B
        this_user_channel = bituser_channel(:,b);
        interference = sigma;
        for i = 1:B
            if i~=b
                interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
            end
        end
        n(b) = sqrt(1+n(b)).*abs(real(this_user_channel'*bituser_precoding(:,b)))/interference;
    end
end
function [gamma,y,z] = update_sinr_term(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T)
    gamma = zeros(T,1);
    y = zeros(B,1);
    z = zeros(B,1);
    power = sum(sum(abs(bituser_precoding).^2))+sum(sum(abs(semuser_precoding).^2));
    %update gamma
    for t = 1:T
        this_user_channel = semuser_channel(:,t);
        interference = power*sigma;
        for i = 1:T
            if i~=t
                interference = interference + abs(this_user_channel'*semuser_precoding(:,i)).^2;
            end
        end
        for i = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        gamma(t) = abs(this_user_channel'*semuser_precoding(:,t)).^2/interference;
    end
    %update y
    for b = 1:B
        this_user_channel = bituser_channel(:,b);
        interference = power*sigma;
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
        interference = power*sigma;
        for i = 1:B
            if i~=b
                interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
            end
        end
        z(b) = abs(this_user_channel'*bituser_precoding(:,b)).^2/interference;
    end
end

function [x,m,n] = update_fraction_variable(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,gamma,y,z,a,b,c,alpha,B,T)
    x = zeros(T,1); 
    m = zeros(B,1);
    n = zeros(B,1);
    coef = alpha*gamma.^(1-alpha);%c*(1-alpha)*gamma.^(alpha)+1;
    power = sum(sum(abs(bituser_precoding).^2))+sum(sum(abs(semuser_precoding).^2));
    %update x
    for t = 1:T
        this_user_channel = semuser_channel(:,t);
        interference = coef(t)*sigma*power+(c+(1-alpha)/(gamma(t))^alpha)*abs(this_user_channel'*semuser_precoding(:,t)).^2;
        for i = 1:T
            if i~=t
                interference = interference + coef(t)*abs(this_user_channel'*semuser_precoding(:,i)).^2;
            end
        end
        for i = 1:B
            interference = interference + coef(t)*abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        x(t) = abs(real(this_user_channel'*semuser_precoding(:,t))) / interference;
    end
    %update m
    for b = 1:B
        this_user_channel = bituser_channel(:,b);
        interference = power*sigma;
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
        interference = power*sigma;
        for i = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        n(b) = sqrt(1+z(b)).*abs(real(this_user_channel'*bituser_precoding(:,b)))/interference;
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
function [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda,mu)
    bituser_matrix = zeros(Nt);
    eye_scalar = 0;
    % inverse matrix calculation
    for bit_idx = 1:B
        scalar = (m(bit_idx)^2*lamda(bit_idx)*J/M + n(bit_idx)^2*lamda(bit_idx)*(L-J))/M;
        eye_scalar = eye_scalar + scalar;
        bituser_matrix = bituser_matrix + scalar* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
    end
    for sem_idx = 1:T
        slope = b*alpha*gamma(sem_idx)^(1-alpha);%b*alpha*gamma(sem_idx)^(alpha-1);
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
        slope = b*alpha*gamma(sem_idx)^(1-alpha);%b*alpha*gamma(sem_idx)^(alpha-1);
        scalar = slope * x(sem_idx)^2;
        %eye_scalar = eye_scalar + scalar;
        semuser_matrix = semuser_matrix + scalar * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    semuser_matrix = semuser_matrix + eye_scalar * sigma * eye(Nt);
    inverse_semuser_matrix = semuser_matrix^(-1);

    % semuser precoding recovery
    reconstruct_semuser_precoding = zeros(Nt,T);
    for sem_idx = 1:T
        %coef = max(c*(1-alpha)*gamma(sem_idx).^(alpha)+1,0.01);
        %slope = b*alpha*gamma(sem_idx)^(alpha-1)/coef;
        coef = max(c+(1-alpha)/gamma(sem_idx)^(alpha),0.01);
        slope = b;
        scalar = slope * x(sem_idx); 
        %scalar2 = (b*alpha*gamma(sem_idx)^(1-alpha)*(1 -c*alpha*gamma(sem_idx)^(alpha-1)/coef)) * x(sem_idx)^2;
        scalar2 = b*(alpha*gamma(sem_idx)^(1-alpha)-coef)* x(sem_idx)^2;
        result_matrix = inverse_semuser_matrix + scalar2/(1-scalar2*semuser_channel(:,sem_idx)'*inverse_semuser_matrix*semuser_channel(:,sem_idx))*inverse_semuser_matrix*semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)'*inverse_semuser_matrix;
        precoding_instance = scalar * result_matrix* semuser_channel(:,sem_idx);
        reconstruct_semuser_precoding(:,sem_idx) = precoding_instance;
    end
end
function [bituser_precoding,semuser_precoding,solving_result,obj] = update_precoding_matrix(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M)
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
        [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,reconstruct_lamda,reconstruct_mu);
        snr_scalar = sum(sum(abs(reconstruct_bituser_precoding).^2))+sum(sum(abs(reconstruct_semuser_precoding).^2));
        diff_factor = (snr_scalar-1)*sigma*(J/M*(m.^2)+(L-J)/M*(n.^2));
        reconstruct_lamda = update_lambda(u+diff_factor,v,L,J,M,m,n,bituser_channel,semuser_channel,reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,B,T);
        iter = iter + 1;
        if iter > max_iter 
            break;
        end
    end
    [bituser_precoding,semuser_precoding,~,~]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,reconstruct_lamda,reconstruct_mu);
    solving_result = 1;
    [obj_list,~,~] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T);
    obj = sum(obj_list);
    reconstruct_lamda
end
function [bituser_precoding,semuser_precoding,solving_result,obj] = update_precoding_matrix_cvx(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M)
    cvx_begin quiet
        variable bituser_precoding(Nt,B) complex;
        variable semuser_precoding(Nt,T) complex;
        dual variable lamda
        dual variable mu        
        expression qos(B,1);
        obj = 0;
        coef = max(c*(1-alpha)*gamma.^(alpha)+1,0.01);
        %calculate the power 
        power = 0;
        for i = 1:T
            power = power + sum_square_abs(semuser_precoding(:,i));
        end
        for i = 1:B
            power = power + sum_square_abs(bituser_precoding(:,i));
        end
        for t = 1:T
            this_user_channel = semuser_channel(:,t);
            interference = coef(t)*sigma*power + c*alpha*gamma(t)^(alpha-1)*pow_abs(this_user_channel'*semuser_precoding(:,t),2);
            slope = b*alpha*gamma(t)^(alpha-1)/(coef(t));
            for i = 1:T
                if i ~= t
                    interference = interference+coef(t)*pow_abs(this_user_channel'*semuser_precoding(:,i),2);
                end
            end
            for i = 1:B
                interference = interference+coef(t)*pow_abs(this_user_channel'*bituser_precoding(:,i),2);
            end
            stationay_point = a+b*(1-alpha)*gamma(t)^(alpha)/coef(t);
%             obj = obj + (stationay_point + slope*(2*x(t)*real(this_user_channel'*semuser_precoding(:,t))-x(t)^2*interference -gamma(t)));
            %obj = min(obj,min(stationay_point + slope*(2*x(t)*real(this_user_channel'*semuser_precoding(:,t))-x(t)^2*interference -gamma(t)),a+b/c));
            obj = obj + min(stationay_point + slope*(2*x(t)*real(this_user_channel'*semuser_precoding(:,t))-x(t)^2*interference),a+b/c);
            %obj = obj + (stationay_point + slope*(2*x(t)*real(this_user_channel'*semuser_precoding(:,t))-x(t)^2*interference));
        end
        

        for b = 1:B
            this_user_channel = bituser_channel(:,b);
            interference_m = power*sigma; interference_n = power*sigma;
            for i = 1:T
                interference_m = interference_m + pow_abs(this_user_channel'*semuser_precoding(:,i),2);
            end
            for i = 1:B
                interference_m = interference_m + pow_abs(this_user_channel'*bituser_precoding(:,i),2);
                interference_n = interference_n + pow_abs(this_user_channel'*bituser_precoding(:,i),2);
            end
            qos(b) = J/M*(log2(1+y(b))-y(b)+2*m(b)*sqrt(1+y(b))*real(this_user_channel'*bituser_precoding(:,b))-m(b)^2*interference_m);
            qos(b) = qos(b) + (L-J)/M*(log2(1+z(b))-z(b)+2*n(b)*sqrt(1+z(b))*real(this_user_channel'*bituser_precoding(:,b))-n(b)^2*interference_n);
            qos(b) = qos(b) - beta(b);
        end
        maximize(obj)
        subject to
           lamda: qos >= 0
           %mu: power <= 1    
    cvx_end
    bituser_precoding = bituser_precoding;
    semuser_precoding = semuser_precoding;
    if strcmp(cvx_status,'Solved')==1  
        solving_result = 1;
    else
        solving_result = 0;
    end
end