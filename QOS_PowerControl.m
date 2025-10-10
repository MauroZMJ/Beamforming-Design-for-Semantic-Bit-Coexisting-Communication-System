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
            [y, z] = update_sinr_term(hv_matrix,p_bit,p_sem,sigma,B,T);
            [x,m, n] = update_fraction_variable(hv_matrix,p_bit,p_sem,sigma,y,z,B,T);
            bituser_precoding = precoding_matrix(:,1:B)*diag(sqrt(p_bit));  semuser_precoding = precoding_matrix(:,B+1:K)*diag(sqrt(p_sem));
            [object_performance,qos,power] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T);
            [p_bit,p_sem,solving_result,obj_value] = update_precoding_matrix_sum(hv_matrix,sigma,Nt,B,T,y,z,x,m,n,beta,L,J,M);
            if abs(mean(object_performance)-last_performance) <= 1e-3 && max(abs(qos-beta)) <= 1e-2
                end_condition = 1;
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
function [y, z] = update_sinr_term(hv_matrix,p_bit,p_sem,sigma,B,T)
    y = zeros(B,1);
    z = zeros(B,1);
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
function [x,m,n] = update_fraction_variable(hv_matrix,p_bit,p_sem,sigma,y,z,B,T)
    x = zeros(T,1);
    m = zeros(B,1);
    n = zeros(B,1);
    %update x
    for t = 1:T
        interference = sigma;
        for i = 1:T
            if i ~=t
                interference = interference + p_sem(i)*abs(hv_matrix(B+t,B+i)).^2; 
            end
        end
        for i = 1:B
            interference = interference + p_bit(i)*abs(hv_matrix(B+t,i)).^2;
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

function [p_bit,p_sem,solving_result,obj_value] = update_precoding_matrix_sum(hv_matrix,sigma,Nt,B,T,y,z,x,m,n,beta,L,J,M)
    cvx_begin quiet
    variable p_bit(B,1);
    variable p_sem(T,1);
    variable o(1);
    %expression obj_list(T,1);
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
    obj_list = 0;
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
        obj_list = obj_list + 2*x(t)*sqrt(p_sem(t))*real(hv_matrix(B+t,B+t))- x(t)^2 * interference ;%- o;
    end
    %maximize(o)
    maximize(sum(obj_list))
    subject to
    power <= 1
    qos_list >= 0
    %obj_list >= 0
    p_sem >= 0
    p_bit >= 0
    cvx_end
    if strcmp(cvx_status,'Solved')==1
        solving_result = 1;
    else
        solving_result = 0;
    end
    obj_value = sum(obj_list);
end