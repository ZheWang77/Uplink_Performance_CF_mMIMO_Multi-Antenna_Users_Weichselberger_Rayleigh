function [HH,R_Vec,H_Vec,Omega,U_APR,U_UET] = functionChannelGeneration_Only_Row(channelGain,M,K,N,L,nbrOfRealizations)
%%=============================================================
%The file is used to generate the coupling matrix with only one row (case (d) in Fig.2) based channel of the paper:
%
% Z. Wang, J. Zhang, B. Ai, C. Yuen and M. Debbah, "Uplink Performance of Cell-Free Massive MIMO With Multi-Antenna Users 
% Over Jointly-Correlated Rayleigh Fading Channels," in IEEE Transactions on Wireless Communications, 
% vol. 21, no. 9, pp. 7391-7406, Sep. 2022, doi: 10.1109/TWC.2022.3158353.

%
%Download article: https://arxiv.org/abs/2110.04962 or https://ieeexplore.ieee.org/document/9737367/
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%%=============================================================

%Prepare to store the results   
U_APR = zeros(M*L,M*L,K);
U_UET = zeros(K*N,K*N,M);

a = L*N/(2*(L*N-1));
% Omega0 = a*ones(L,N);
% Omega0(1,1) = L*N/2;

Omega = zeros(M*L,K*N);
WW = sqrt(0.5)*(randn(M*L,K*N,nbrOfRealizations)+1i*randn(M*L,K*N,nbrOfRealizations));
R_Vec = zeros(L*N,L*N,M,K);
H_Vec = zeros(L*N,nbrOfRealizations,M,K);
HH = zeros(M*L,nbrOfRealizations,K*N);


for k = 1:K
    for m = 1:M
        
        %----Unitary Matrix
        [U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k),~,~] = svd(rand(L,L) + 1i*rand(L,L));
        [U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m),~,~] = svd(rand(N,N) + 1i*rand(N,N)); 

        %----Coupling Matrix
        
        Omega0 = sort(rand(N,1),'descend');
        Omega((m-1)*L+1,(k-1)*N+1:k*N) = channelGain(m,k)*Omega0/sum(Omega0(:))*L*N; % generate the coupling matrix randomly
        
           
    end
end


%---Channel Generation
H = permute((sqrt(Omega).*WW),[1,3,2]);
for k = 1:K
    for i = (k-1)*N+1:k*N
        
        HH(:,:,i) = U_APR(:,:,k)*H(:,:,i);
        
    end
end

for m = 1:M
    for t = (m-1)*L+1:m*L
        
        H11 = reshape(HH(t,:,:),nbrOfRealizations,K*N);
        HH(t,:,:) = reshape(H11*U_UET(:,:,m)',nbrOfRealizations,K*N);
        
    end
end



%---Full Correlation Matrix & Channel Vectorization
HH_v = permute(HH,[1,3,2]);
for m = 1:M
    for k = 1:K
        
        Omeg = reshape(Omega((m-1)*L+1:m*L,(k-1)*N+1:k*N),L*N,1);
        R_Vec(:,:,m,k) = kron(conj(U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m)),U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k))*diag(Omeg)*kron(conj(U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m)),U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k))';
        H_Vec(:,:,m,k) = reshape(HH_v((m-1)*L+1:m*L,(k-1)*N+1:k*N,:),L*N,nbrOfRealizations);
        
    end
end
