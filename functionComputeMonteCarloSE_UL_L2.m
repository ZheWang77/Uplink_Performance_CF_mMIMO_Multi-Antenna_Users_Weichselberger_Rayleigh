function [SE] = functionComputeMonteCarloSE_UL_L2(V,H,tau_c,tau_p,nbrOfRealizations,M,K,L,N,p)
%%=============================================================
%This function is used to compute the achievable SE for the Level 2 of the paper:
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

%Compute the pre-log factor
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store simulation results
SE = zeros(K,1);
r = zeros(N,N,K);
pn = kron(p,ones(N,1));
Term1x = zeros(N,N,K,nbrOfRealizations);
Term2x = zeros(N,N,K,nbrOfRealizations);
Term3x = zeros(N,N,K,nbrOfRealizations);


% Go through all channel realizations
for i = 1:nbrOfRealizations
    
    %All channels for one realization
    Hallj = reshape(H(:,i,:),[M*L K*N]);
    Vallj = reshape(V(:,i,:),[M*L K*N]);
    
    for k = 1:K
        
        Term1x(:,:,k,i) = Vallj(:,(k-1)*N+1:k*N)'*Hallj(:,(k-1)*N+1:k*N);
        Term2x(:,:,k,i) = diag(pn((k-1)*N+1:k*N))*(Vallj(:,(k-1)*N+1:k*N)'*Hallj)*(Vallj(:,(k-1)*N+1:k*N)'*Hallj)';
        Term3x(:,:,k,i) = Vallj(:,(k-1)*N+1:k*N)'*Vallj(:,(k-1)*N+1:k*N);
        
    end
end


Term1 = mean(Term1x,4);
Term2 = mean(Term2x,4);
Term3 = mean(Term3x,4);

for k = 1:K
    
    r(:,:,k) = eye(N) + p(k)*Term1(:,:,k)'/((Term2(:,:,k)) - p(k)*Term1(:,:,k)*Term1(:,:,k)' + Term3(:,:,k))*Term1(:,:,k);

end
    

%Calculate the SE of each UE k
for k = 1:K
    
    SE(k) = prelogFactor*real(log2(det(r(:,:,k))));

end

