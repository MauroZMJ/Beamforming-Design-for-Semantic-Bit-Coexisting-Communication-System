function [semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = cal_SINR(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T)
    semuser_sinr = zeros(T,1);
    bituser_bit_sinr = zeros(B,1);
    bituser_sem_sinr = zeros(B,1);
    for t = 1:T
        interference = sigma;
        this_user_channel = semuser_channel(:,t);
        for i = 1:T
            if i~=t
                interference = interference + abs(this_user_channel'*semuser_precoding(:,i)).^2;
            end
        end
        for i = 1:B
            interference = interference + abs(this_user_channel'*bituser_precoding(:,i)).^2;
        end
        semuser_sinr(t) = abs(this_user_channel'*semuser_precoding(:,t)).^2/interference;
    end
    for b = 1:B
        interference_a = sigma; 
        this_user_channel = bituser_channel(:,b);
        for i = 1:B
            if i~=b
                interference_a = interference_a + abs(this_user_channel'*bituser_precoding(:,i)).^2;
            end
        end
        interference_b = interference_a;
        for i = 1:T
            interference_b = interference_b + abs(this_user_channel'*semuser_precoding(:,i)).^2;
        end
        bituser_bit_sinr(b) = abs(this_user_channel'*bituser_precoding(:,b)).^2/interference_a;
        bituser_sem_sinr(b) = abs(this_user_channel'*bituser_precoding(:,b)).^2/interference_b;       
    end
end