clear all
load('workspace.mat')
b = 2.37; %mu=0.6544;

%[reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda,sigma);
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
        snr_scalar = sum(sum(abs(reconstruct_bituser_precoding).^2))+sum(sum(abs(reconstruct_semuser_precoding).^2));
        diff_factor = (snr_scalar-1)*sigma*(J/M*(m.^2)+(L-J)/M*(n.^2));
        reconstruct_lamda = update_lambda(u+diff_factor,v,L,J,M,m,n,bituser_channel,semuser_channel,reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,B,T);
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
        slope = b*alpha*gamma(sem_idx)^(alpha-1);
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
        slope = b*alpha*gamma(sem_idx)^(alpha-1);
        scalar = slope * x(sem_idx)^2;
        %eye_scalar = eye_scalar + scalar;
        semuser_matrix = semuser_matrix + scalar * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
    end
    semuser_matrix = semuser_matrix + eye_scalar * sigma * eye(Nt);
    inverse_semuser_matrix = semuser_matrix^(-1);

    % semuser precoding recovery
    reconstruct_semuser_precoding = zeros(Nt,T);
    for sem_idx = 1:T
        slope = b*alpha*gamma(sem_idx)^(alpha-1)/(1+c*(1-alpha)*gamma(sem_idx)^(alpha));
        scalar = slope * x(sem_idx); 
        scalar2 = (b*alpha*gamma(sem_idx)^(alpha-1)*(1 -c*alpha*gamma(sem_idx)^(alpha-1)/(1+(c*(1-alpha)*gamma(sem_idx)^(alpha))))) * x(sem_idx)^2;
        result_matrix = inverse_semuser_matrix + scalar2/(1-scalar2*semuser_channel(:,sem_idx)'*inverse_semuser_matrix*semuser_channel(:,sem_idx))*inverse_semuser_matrix*semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)'*inverse_semuser_matrix;
        precoding_instance = scalar * result_matrix* semuser_channel(:,sem_idx);
        reconstruct_semuser_precoding(:,sem_idx) = precoding_instance;
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