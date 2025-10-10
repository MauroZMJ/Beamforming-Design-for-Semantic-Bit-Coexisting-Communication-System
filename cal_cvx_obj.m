function obj = cal_cvx_obj(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M)
    obj = zeros(T,1);
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
        obj(t) = stationay_point + slope*(2*x(t)*real(this_user_channel'*semuser_precoding(:,t))-x(t)^2*interference -gamma(t));
    end
end