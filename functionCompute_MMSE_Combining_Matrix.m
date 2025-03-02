function [V_MMSE_Combining] = functionCompute_MMSE_Combining_Matrix(Hhat,C,nbrOfRealizations,L,N,K,M,p)
%%=============================================================
%This function is used to generate the L-MMSE combining matrix of the paper:
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



%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
   p = p*ones(K,1);
end


pn = kron(p,ones(1,N));
%Store identity matrices of different sizes
eyeL = eye(L);
V_MMSE_Combining = zeros(M*L,nbrOfRealizations,K*N);


%Compute sum of all estimation error correlation matrices at every BS
C_tot = sum(C,4);


%Diagonal matrix with transmit powers and its square root
Dp = diag(pn);

%% Go through all channel realizations
for n = 1:nbrOfRealizations

    %Go through all APs
    for m = 1:M
        
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(m-1)*L:m*L,n,:),[L K*N]);

        %Compute MMSE combining
        V_MMSE_Combining((m-1)*L+1:m*L,n,:) = ((Hhatallj*Dp*Hhatallj') + C_tot(:,:,m)+eyeL)\(Hhatallj*Dp);
        
    end
end