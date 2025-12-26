function [object_performance,qos,power] = performance(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,a,b,c,alpha,L,J,M,B,T)
    power = sum(sum(abs(bituser_precoding).^2))+sum(sum(abs(semuser_precoding).^2));
    [semuser_sinr,bituser_bit_sinr,bituser_sem_sinr] = cal_SINR(bituser_channel,semuser_channel,bituser_precoding,semuser_precoding,sigma,B,T);
    qos = J/M*(log2(1+bituser_sem_sinr)) + (L-J)/M*(log2(1+bituser_bit_sinr));
    object_performance = a + b./(c+semuser_sinr.^(-alpha));
end