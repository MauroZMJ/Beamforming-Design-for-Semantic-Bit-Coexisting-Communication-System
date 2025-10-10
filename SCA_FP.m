function [precoding_matrix,object_performance_list,cvx_obj_list,solving_result] = SCA_FP(channel,precoding_matrix,B,T,L,J,M,a,b,c,alpha,beta,sigma)
Nt = size(channel,1); K = size(channel,2);
bituser_channel = channel(:,1:B); semuser_channel = channel(:,B+1:K); %seperate the channel of bituser and semuser
bituser_precoding = precoding_matrix(:,1:B); semuser_precoding = precoding_matrix(:,B+1:K);
end_condition = 0;
iteration = 0;
object_performance_list = [];
cvx_obj_list = [];
last_performance = 0;
while end_condition == 0
    [gamma, y, z] = update_sinr_term(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T);
    [x, m, n] = update_fraction_variable(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,y,z,B,T);
    [object_performance,qos,power] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T)
    %cvx_obj = cal_cvx_obj(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M);
    [bituser_precoding,semuser_precoding,solving_result,cvx_obj] = update_precoding_matrix(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M);
    %[semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = cal_SINR(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T)
    iteration = iteration + 1;
    object_performance_list = [object_performance_list,object_performance];
    cvx_obj_list = [cvx_obj_list,cvx_obj];
    if iteration == 10
        la = 1;
        %plot(sum(object_performance_list,1)); hold on
        %plot(sum(cvx_obj_list,1))
    end
end
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

function [x,m,n] = update_fraction_variable(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,y,z,B,T)
    x = zeros(T,1); 
    m = zeros(B,1);
    n = zeros(B,1);
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

function [bituser_precoding,semuser_precoding,solving_result,obj] = update_precoding_matrix(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M)
    cvx_begin
        variable bituser_precoding(Nt,B) complex;
        variable semuser_precoding(Nt,T) complex;
        expression qos(B,1);
        obj = 0;
        for t = 1:T
            interference = sigma;
            this_user_channel = semuser_channel(:,t);
            slope = alpha*b*gamma(t)^(-alpha-1)/(c+gamma(t)^(-alpha))^2;
            for i = 1:T
                if i ~= t
                    interference = interference+pow_abs(this_user_channel'*semuser_precoding(:,i),2);
                end
            end
            for i = 1:B
                interference = interference+pow_abs(this_user_channel'*bituser_precoding(:,i),2);
            end
            stationay_point = a + b/(c+gamma(t)^(-alpha));
%             obj = obj + (stationay_point + slope*(2*x(t)*real(this_user_channel'*semuser_precoding(:,t))-x(t)^2*interference -gamma(t)));
            %obj = min(obj,min(stationay_point + slope*(2*x(t)*real(this_user_channel'*semuser_precoding(:,t))-x(t)^2*interference -gamma(t)),a+b/c));
            obj = obj + min(stationay_point + slope*(2*x(t)*real(this_user_channel'*semuser_precoding(:,t))-x(t)^2*interference -gamma(t)),a+b/c);
        end
        
        %calculate the power 
        power = 0;
        for i = 1:T
            power = power + sum_square_abs(semuser_precoding(:,i));
        end
        for i = 1:B
            power = power + sum_square_abs(bituser_precoding(:,i));
        end
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
            qos(b) = J/M*(log2(1+y(b))-y(b)+2*m(b)*sqrt(1+y(b))*real(this_user_channel'*bituser_precoding(:,b))-m(b)^2*interference_m);
            qos(b) = qos(b) + (L-J)/M*(log2(1+z(b))-z(b)+2*n(b)*sqrt(1+z(b))*real(this_user_channel'*bituser_precoding(:,b))-n(b)^2*interference_n);
            qos(b) = qos(b) - beta(b);
        end
        maximize(obj)
        subject to
           qos >= 0
           power <= 1    
    cvx_end
    bituser_precoding = bituser_precoding;
    semuser_precoding = semuser_precoding;
    if strcmp(cvx_status,'Solved')==1  
        solving_result = 1;
    else
        solving_result = 0;
    end
end