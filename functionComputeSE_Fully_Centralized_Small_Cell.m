function [SE_level4,SE_l1] = functionComputeSE_Fully_Centralized_Small_Cell(H_Combining,H_hat,C,tau_c,tau_p,nbrOfRealizations,N,L,K,M,p,V_d)
warning('off');
%%=============================================================
%This function is used to compute the achievable SE for Level 4 and Level 1 of the paper:
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

%Diagonal matrix with transmit powers and its square root
% Dp12 = diag(sqrt(pn));
Dp = diag(pn);

%Store identity matrices of different sizes
eyeL = eye(L);
eyeML = eye(M*L);

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
% SE_l1 = zeros(K);

%Compute sum of all estimation error correlation matrices at every BS
C_tot = zeros(L,L,M);

for k = 1:K
    C_tot = C_tot + p(k)*C(:,:,:,k);
end

C_tot_blk = zeros(M*L,M*L);

for m = 1:M
    
    C_tot_blk(1+(m-1)*L:m*L,1+(m-1)*L:m*L) = C_tot(:,:,m);
    
end


%Prepare to save simulation results
SE_level1 = zeros(K,M);
SE_level4 = zeros(K,1);
%% Go through all channel realizations
for i = 1:nbrOfRealizations
    
    
     %--------------Level 4
    
    %Extract channel estimate realizations from all UEs to all APs
    Hhatallj = reshape(H_hat(:,i,:),[M*L K*N]);
    
    if V_d == 1
        
        %Compute MR combining
        V = Hhatallj;
        
    else
        %Compute C-MMSE combining
        V = ((Hhatallj*Dp*Hhatallj')+C_tot_blk+eyeML)\(Hhatallj*Dp);
    end
        
 
    %Go through all UEs
    for k = 1:K
        
        v = V(:,(k-1)*N+1:k*N); %Extract combining matrix
        
         %Compute numerator and denominator of instantaneous SINR at Level 4
        numerator = sqrt(p(k))*(v'*Hhatallj(:,(k-1)*N+1:k*N));
        denominator = v'*(Hhatallj*Dp*Hhatallj' + C_tot_blk + eyeML)*v - numerator*numerator';
        
        SE_level4(k) = SE_level4(k) + prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)))/nbrOfRealizations;
            

    end

       
  %----Level1
    
    %Go through all APs
    for m = 1:M
        
        
        %Extract channel estimate realizations from all UEs to AP m
        Hhatallj1 = reshape(H_Combining(1+(m-1)*L:m*L,i,:),[L K*N]);
        Hallj1 = reshape(H_hat(1+(m-1)*L:m*L,i,:),[L K*N]);

        %Compute MR/L-MMSE combining
        V = Hhatallj1;

        %Go through all UEs
        for k = 1:K
            
            v = V(:,(k-1)*N+1:k*N); %Extract combining matrix
            numerator = sqrt(p(k))*v'*Hallj1(:,(k-1)*N+1:k*N);
            denominator = v'*(Hallj1*Dp*Hallj1' + C_tot(:,:,m) + eyeL)*v - numerator*numerator';
            
            SE_level1(k,m) = SE_level1(k,m) + prelogFactor*real(log2(det(eye(N)+numerator'/denominator*numerator)))/nbrOfRealizations;
        end
    end
end


%Compute SE for Level 1
SE_l1 = max(SE_level1,[],2);
