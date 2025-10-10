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