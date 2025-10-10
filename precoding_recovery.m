function [reconstruct_bituser_precoding,reconstruct_semuser_precoding,inverse_bituser_matrix,inverse_semuser_matrix]=precoding_recovery(bituser_channel,semuser_channel,sigma,Nt,B,T,gamma,y,z,x,m,n,beta,alpha,a,b,c,L,J,M,lamda,mu)
bituser_matrix = eye(Nt)*mu;
%% inverse matrix calculation
for bit_idx = 1:B
    scalar = (m(bit_idx)^2*lamda(bit_idx)*J/M + n(bit_idx)^2*lamda(bit_idx)*(L-J))/M;
    bituser_matrix = bituser_matrix + scalar* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
end
for sem_idx = 1:T
    slope = b*alpha*gamma(sem_idx)^(alpha-1);
    scalar = slope * x(sem_idx)^2;
    bituser_matrix = bituser_matrix + scalar * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
end
inverse_bituser_matrix = bituser_matrix^(-1);
%% bituser precoding recovery
reconstruct_bituser_precoding = zeros(Nt,B);
for bit_idx = 1:B
    scalar = lamda(bit_idx)*m(bit_idx)*sqrt(1+y(bit_idx))*J/M + lamda(bit_idx)*n(bit_idx)*sqrt(1+z(bit_idx))*(L-J)/M;
    precoding_instance = scalar * inverse_bituser_matrix * bituser_channel(:,bit_idx);
    reconstruct_bituser_precoding(:,bit_idx) = precoding_instance;
end

%% inverse matrix calculation
semuser_matrix = eye(Nt)*mu;
for bit_idx = 1:B
    scalar = lamda(bit_idx)*J/M*m(bit_idx)^2;
    semuser_matrix = semuser_matrix + scalar* bituser_channel(:,bit_idx)*bituser_channel(:,bit_idx)';
end
for sem_idx = 1:T
    slope = b*alpha*gamma(sem_idx)^(alpha-1);
    scalar = slope * x(sem_idx)^2;
    semuser_matrix = semuser_matrix + scalar * semuser_channel(:,sem_idx)*semuser_channel(:,sem_idx)';
end
inverse_semuser_matrix = semuser_matrix^(-1);

%% semuser precoding recovery
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