clear all
load('workspace.mat')
b = 2.37; %mu=0.6544;

%[reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda,mu);
u = zeros(B,1); v = zeros(B,1); reconstruct_lamda = ones(B,1)*0.01; reconstruct_lamda_old = zeros(B,1);
for bit_idx = 1:B
    u(bit_idx) = beta(bit_idx) - (J/M)*(log2(1+y(bit_idx))-y(bit_idx)-m(bit_idx)^2*sigma) - ((L-J)/M)*(log2(1+z(bit_idx))-z(bit_idx)-n(bit_idx)^2*sigma);
    v(bit_idx) = 2*(J/M*m(bit_idx)*sqrt(1+y(bit_idx))+(L-J)/M*n(bit_idx)*sqrt(1+z(bit_idx)))^2;
end
i = 0
reconstruct_mu = 1; reconstruct_mu_old = 0; 
tic
while sum(abs(reconstruct_mu-reconstruct_mu_old)) >1e-6
    %fixed algorithm for the search of lamda
    j=0;
    while sum(abs(reconstruct_lamda-reconstruct_lamda_old))>1e-5
        reconstruct_lamda_old = reconstruct_lamda;
        [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,reconstruct_lamda,reconstruct_mu);
        reconstruct_lamda = update_lambda(u,v,L,J,M,m,n,bituser_channel,semuser_channel,reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,B,T);
        %sum(abs(reconstruct_lamda-reconstruct_lamda_old))
        j = j + 1;
    end
    reconstruct_lamda_old = zeros(B,1);
    reconstruct_mu_old = reconstruct_mu;
    reconstruct_mu = update_mu(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,reconstruct_lamda)
    i = i +1;
    j
end
toc



