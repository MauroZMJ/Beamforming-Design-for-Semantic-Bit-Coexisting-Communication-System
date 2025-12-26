function [H_list] = generate_channel(data_num,Nt,K,Lp)
% data_num: number of channel data to generate
% Nt: number of transmit antennas
% K: number of users
% Lp: number of path
% Return H: Downlink channel, Nt x K, complex 
    H_list = zeros(data_num,Nt,K)+1j*zeros(data_num,Nt,K); 
    for num=1:1:data_num
        if mod(num,1000)==0
            sprintf('%d data generated',num)
        end
        H = zeros(Nt,K)+1j*zeros(Nt,K);
        for k=1:1:K
            AoDs = rand(Lp,1)*2*pi; % angle of departure at BS, uniform distribution in [0,2*pi]
            Gains = sqrt(1/2)*(randn(Lp,1)+1j*randn(Lp,1));% channel gains of paths
            Hk = zeros(Nt,1)+1j*zeros(Nt,1);
            a_t = zeros(Nt,1)+1j*zeros(Nt,1);
            for lp=1:1:Lp
                for nt=1:1:Nt
                    a_t(nt) = exp(-1j*2*pi*1/2*(nt-1)*sin(AoDs(lp))); 
                end
                Hk = Hk + Gains(lp)*a_t;
            end
            H(:,k) = Hk/sqrt(Lp);
        end
        H_list(num,:,:) = H;
    end
end

