function [Phi,Omega,C_MMSE,C_MMSE_MMSE_Combining] = functionMatrixGeneration(R_Vec,M,K,L,N,tau_p,p,Pset)
%%=============================================================
%This function is used to generate the matrix used in the next
%subsequent calculation of the paper:
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


%Prepare to store the result
Phi = zeros(L*N,L*N,M,K);
Omega = zeros(L*N,L*N,M,K);
C_MMSE = zeros(L*N,L*N,M,K);
C_MMSE_MMSE_Combining = zeros(L,L,M,K); %Matrix for MMSE Combining
 
%Go through all APs          
for m = 1:M
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the UEs indexes that use the same pilot as UE k
        inds = Pset(:,k);
        PsiInv = zeros(L*N,L*N);
        
        %Go through all UEs that use the same pilot as UE k 
        for z = 1:length(inds)   
            
            PsiInv = PsiInv + p(inds(z))*tau_p*R_Vec(:,:,m,inds(z));


        end
            PsiInv = PsiInv + eye(L*N);

            
            for z = 1:length(inds)
                
                Phi(:,:,m,inds(z)) = PsiInv;
            
            end
            
            Omega(:,:,m,k) = p(k)*tau_p*R_Vec(:,:,m,k)/PsiInv*R_Vec(:,:,m,k);
            
    end
end

%Generate estimation error correlation matrices
for k = 1:K
    
    C_MMSE(:,:,:,k) = R_Vec(:,:,:,k) - Omega(:,:,:,k);
   
end

for k = 1:K
    for m = 1:M
        for n = 1:N
            
            C_MMSE_MMSE_Combining(:,:,m,k) = C_MMSE_MMSE_Combining(:,:,m,k) + C_MMSE((n-1)*L+1:n*L,(n-1)*L+1:n*L,m,k);
            
        end
    end
end
        
            
