function V = ZF_beamforming(channel,B,T,L,J,M,a,b,c,alpha,beta,sigma)
    Nt = size(channel,1); K = B+T;
    V = zeros(Nt,B+T);
    H_instance = channel(:,1:K);
    V_instance = H_instance*(H_instance'*H_instance)^(-1);
    for i = 1:K
       V(:,i) = V_instance(:,i)/norm(V_instance(:,i),2)*exp(-1i*angle(H_instance(:,i)'*V_instance(:,i))); 
    end
    V = QOS_PowerControl_v4(channel,V,B,T,L,J,M,a,b,c,alpha,beta,sigma);
end  