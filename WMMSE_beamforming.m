function [V] = WMMSE_beamforming(channel_miso,B,T,L,J,M,a,b,c,alpha,beta,sigma)
    K = B+T;
    alpha_wmmse = ones(1,K);
    Nt = size(channel_miso,1); 
    Nr = 1; Rb = 1; dk = 1; Power_constrant = 1;
    H = zeros(Nr,Nt,K,Rb);
    H(1,:,:,1) = conj(channel_miso(:,1:K)); 
    V_instance = zeros(Nt,B+T);
    %generate WMMSE argument 
    V_init = channel_miso*(channel_miso'*channel_miso)^(-1); V_init = V_init*(Power_constrant/sum(sum(sum(abs(V_init).^2))))^0.5;
    V = zeros(Nt,dk,K,Rb); V(:,1,:,1) = V_init(:,1:K);
    U = zeros(Nr,dk,K,Rb);
    W = ones(dk,dk,K,Rb);

    new_V = zeros(Nt,dk,K);
    new_U = zeros(Nr,dk,K,Rb);
    new_W = zeros(dk,dk,K,Rb);

    %WMMSE
    object_value=10;
    episilon=1e-5; % 1e-4
    max_iters = 150; % another kind of stop criterion
    iter = 0;
    while object_value>episilon
        iter = iter + 1;
        if iter>=max_iters
            break;
        end
        for rb = 1:Rb
            for user_index=1:K
                %update Uk
                inverse = zeros(Nr,Nr);
                for i = 1:K
                    inverse = inverse+H(:,:,user_index,rb)*V(:,:,i)*V(:,:,i)'*H(:,:,user_index,rb)';
                end
                new_U(:,:,user_index,rb) = (inverse+sigma/Power_constrant*sum(sum(sum(abs(V).^2)))*eye(Nr))^(-1)*H(:,:,user_index,rb)*V(:,:,user_index);
            end
        end
        for rb = 1:Rb
            for user_index=1:K
                %update Wk
                new_W(:,:,user_index,rb) = (eye(dk)-new_U(:,:,user_index,rb)'*H(:,:,user_index,rb)*V(:,:,user_index))^(-1);
            end
        end

        sum_trace_UWU = 0;
        alpha_sum_HUWUH = zeros(Nt,Nt);
        for rb = 1:Rb
            for j=1:K
                sum_trace_UWU = sum_trace_UWU+trace(new_U(:,:,j,rb)*new_W(:,:,j,rb)*new_U(:,:,j,rb)');
                alpha_sum_HUWUH = alpha_sum_HUWUH+alpha_wmmse(j)*H(:,:,j,rb)'*new_U(:,:,j,rb)*new_W(:,:,j,rb)*new_U(:,:,j,rb)'*H(:,:,j,rb);
            end
        end
        for user_index=1:K
            HUW_this_user = zeros(Nt,dk);
            for rb = 1:Rb
                HUW_this_user = HUW_this_user + H(:,:,user_index,rb)'*new_U(:,:,user_index,rb)*new_W(:,:,user_index,rb);
            end
            new_V(:,:,user_index) = alpha_wmmse(user_index)*(alpha_sum_HUWUH+sigma/Power_constrant*sum_trace_UWU*eye(Nt))^(-1)*HUW_this_user; 
        end    
        %calculate the object value 
        for l=1:K
          object_value_list_old(l) = log(det(W(:,:,l)));
          object_value_list_new(l) = log(det(new_W(:,:,l)));
        end
        object_value = abs(sum(object_value_list_old)-sum(object_value_list_new));
        U=new_U;
        V=new_V;
        W=new_W;
    end
    V = V*(Power_constrant/sum(sum(sum(abs(V).^2))))^0.5;
    V_instance(:,1:K) = V(:,1,:,1); V_instance = squeeze(V);
    for i = 1:K
        V_instance(:,i) = V_instance(:,i)/norm(V_instance(:,i),2)*exp(-1i*angle(channel_miso(:,i)'*V_instance(:,i)));
    end
    V = QOS_PowerControl_v4(channel_miso,V_instance,B,T,L,J,M,a,b,c,alpha,beta,sigma);
    
end