function [p_bit,p_sem] = QOS_PowerControl_feasible_point(channel,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma)
    Nt = size(channel,1); K = size(channel,2);
    hv_matrix = channel'*precoding_matrix;
    p_bit = ones(B,1)/(B); p_sem = ones(T,1)*0.001;
    end_condition = 0;
    iteration = 0;
    object_performance_list = [];
    cvx_obj_list = [];
    last_performance = 0;
    bituser_channel = channel(:,1:B); semuser_channel = channel(:,B+1:K); %seperate the channel of bituser and semuser
    while end_condition == 0
        [y, z] = update_sinr_term(hv_matrix,p_bit,p_sem,sigma,B,T);
        [m, n] = update_fraction_variable(hv_matrix,p_bit,p_sem,sigma,y,z,B,T);
        bituser_precoding = precoding_matrix(:,1:B)*diag(sqrt(p_bit));  semuser_precoding = precoding_matrix(:,B+1:K)*diag(sqrt(p_sem));
        [~,qos,~] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T);
        %[bit_qos,sem_sinr] = only_for_test(hv_matrix,p_bit,p_sem,sigma,Nt,B,T,y,z,x,m,n,beta,L,J,M);
        [p_bit,p_sem,~,~] = update_precoding_matrix(hv_matrix,p_sem,sigma,Nt,B,T,y,z,m,n,beta,L,J,M);
        if abs(mean(qos)-last_performance) <= 1e-4 || sum(qos >= beta) == B
            end_condition = 1;
        else
            last_performance = mean(qos);
        end
        iteration = iteration + 1;
    end
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
function [m,n] = update_fraction_variable(hv_matrix,p_bit,p_sem,sigma,y,z,B,T)
    m = zeros(B,1);
    n = zeros(B,1);
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
function [p_bit,p_sem,solving_result,obj_value] = update_precoding_matrix(hv_matrix,p_sem,sigma,Nt,B,T,y,z,m,n,beta,L,J,M)
    cvx_begin quiet
    variable p_bit(B,1);
    variable o(1);
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
        qos_list(b) = qos_list(b)* beta(b) - o;
    end

    %calculate the power
    power = 0;
%     for i = 1:T
%         power = power + p_sem(i);
%     end
    for i = 1:B
        power = power + p_bit(i);
    end
    maximize(o)
    subject to
    power <= 1
    qos_list >= 0
    p_bit >= 0
    cvx_end
    if strcmp(cvx_status,'Solved')==1
        solving_result = 1;
    else
        solving_result = 0;
    end
    obj_value = o;
end