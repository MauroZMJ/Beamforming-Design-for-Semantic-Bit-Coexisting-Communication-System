function mu = update_mu(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda)
    mu_low = 0; mu_high = 1;
    [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda,mu_high);
    while sum(sum(abs(reconstruct_semuser_precoding).^2))+sum(sum(abs(reconstruct_bituser_precoding).^2))>1
        mu_high = mu_high * 2;
        [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda,mu_high);
    end
    while abs(mu_high-mu_low) >1e-4
        mu_avg = (mu_high + mu_low)/2;
        [reconstruct_bituser_precoding,reconstruct_semuser_precoding,~,~]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda,mu_avg);
        if sum(sum(abs(reconstruct_semuser_precoding).^2))+sum(sum(abs(reconstruct_bituser_precoding).^2)) < 1
            mu_high = mu_avg;
        else
            mu_low = mu_avg;
        end
    end
    mu = mu_avg;
end