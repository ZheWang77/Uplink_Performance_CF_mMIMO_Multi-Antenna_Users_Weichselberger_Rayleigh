function [H_Vec,R_Vec] = functionChannelVectorization(H,R_AP,R_UE,M,K,N,L,nbrOfRealizations)
%%=============================================================
%The file is used to implement the vectorization step to the channel matrix
%of the paper:
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



H_Vec = zeros(L*N,nbrOfRealizations,M,K);
R_Vec = zeros(L*N,L*N,M,K);

HH = permute(H,[1,3,2]);

for m = 1:M
    for k = 1:K
        
       H_Vec(:,:,m,k) = reshape(HH((m-1)*L+1:m*L,(k-1)*N+1:k*N,:),L*N,nbrOfRealizations);
       R_Vec(:,:,m,k) = kron(R_UE(:,:,m,k).',R_AP(:,:,m,k));
%        R_Vec(:,:,m,k) = kron(R_UE(:,:,m,k),R_AP(:,:,m,k));
    end
end