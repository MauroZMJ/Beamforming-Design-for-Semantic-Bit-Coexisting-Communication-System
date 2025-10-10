function [channel_dataset,precoding_matrix_dataset] = generate_dataset(data_number,Nt,B,T,Lp,L,J,M,a,b,c,alpha,beta,sigma)
    channel_dataset = generate_channel(data_number,Nt,B+T,Lp);
    precoding_matrix_dataset = zeros(size(channel_dataset));
    for i = 1:data_number
        precoding_matrix_dataset(i,:,:) = SINR_balance_beamforming(squeeze(channel_dataset(i,:,:)),B,T,L,J,M,a,b,c,alpha,beta,sigma,1e-3);
        if mod(i,20) == 0
            i
        end
    end
end