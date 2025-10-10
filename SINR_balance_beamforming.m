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
function [bituser_precoding,semuser_precoding,solving_result,obj_value] = update_precoding_matrix(bituser_channel,semuser_channel,semuser_precoding,sigma,Nt,B,T,y,z,m,n,beta,L,J,M)
    cvx_begin quiet
        variable bituser_precoding(Nt,B) complex;
        variable o(1);
        dual variable lamda
        dual variable nu
        expression obj_list(B,1);
        for b = 1:B
            this_user_channel = bituser_channel(:,b);
            interference_m = sigma; interference_n = sigma;
            for i = 1:T
                interference_m = interference_m + pow_abs(this_user_channel'*semuser_precoding(:,i),2);
            end
            for i = 1:B
                interference_m = interference_m + pow_abs(this_user_channel'*bituser_precoding(:,i),2);
                interference_n = interference_n + pow_abs(this_user_channel'*bituser_precoding(:,i),2);
            end
            obj_list(b) = J/M*(log2(1+y(b))-y(b)+2*m(b)*sqrt(1+y(b))*real(this_user_channel'*bituser_precoding(:,b))-m(b)^2*interference_m);
            obj_list(b) = obj_list(b) + (L-J)/M*(log2(1+z(b))-z(b)+2*n(b)*sqrt(1+z(b))*real(this_user_channel'*bituser_precoding(:,b))-n(b)^2*interference_n);
            obj_list(b) = obj_list(b) - o*beta(b);
        end
        
        %calculate the power
        power = 0;
        for i = 1:T
            power = power + sum_square_abs(semuser_precoding(:,i));
        end
        for i = 1:B
            power = power + sum_square_abs(bituser_precoding(:,i));
        end
        maximize(o)
        subject to
           nu: power <= 1
           lamda: obj_list >= 0
    cvx_end
    bituser_precoding = bituser_precoding;
    semuser_precoding = semuser_precoding;
    if strcmp(cvx_status,'Solved')==1  
        solving_result = 1;
    else
        solving_result = 0;
    end
    obj_value = o;
end


