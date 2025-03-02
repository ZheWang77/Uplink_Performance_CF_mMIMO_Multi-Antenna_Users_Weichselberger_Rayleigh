function [SE] = functionComputeMonteCarloSE_UL_LSFD(H_Combining,H,tau_c,tau_p,nbrOfRealizations,N,L,K,M,p)
%%=============================================================
%This function is used to compute the achievable SE for the LSFD scheme of the paper:
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


%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE = zeros(K,1);
signal_MR = zeros(M*N,N,K);
scaling_MR = zeros(M*N,M*N,K);
Gp_MR_total = zeros(M*N,M*N,K,K);
A = zeros(M*N,N,K);


%Go through all channel realizations
for i = 1:nbrOfRealizations
    
    
    %-----------------Levels 1-3
    gp_MR = zeros(M*N,N,K,K);
    
    %Go through all APs
    for m = 1:M
        
        %Extract channel realizations from all UEs to AP m
        Hallj = reshape(H(1+(m-1)*L:m*L,i,:),[L K*N]);
        %Extract channel estimate realizations from all UEs to AP m
        Hhatallj = reshape(H_Combining(1+(m-1)*L:m*L,i,:),[L K*N]);

        
        %Compute MR/L-MMSE combining
        V = Hhatallj;
        
        for k = 1:K
            
            v = V(:,(k-1)*N+1:k*N); %Extract combining matrix
            signal_MR((m-1)*N+1:m*N,:,k) = signal_MR((m-1)*N+1:m*N,:,k) + (v'*Hallj(:,(k-1)*N+1:k*N))/nbrOfRealizations; % (v_mk)'h_mk
            scaling_MR((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = scaling_MR((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) + v'*v/nbrOfRealizations; % ||v_mk||^2
   
            for l = 1:K
                
                gp_MR((m-1)*N+1:m*N,:,k,l) = v'*Hallj(:,(l-1)*N+1:l*N);
            
            end
        end
    end
    
    
    
    
    for k = 1:K
        for l = 1:K
            

            Gp_MR_total(:,:,k,l) = Gp_MR_total(:,:,k,l) + p(l)*(gp_MR(:,:,k,l)*gp_MR(:,:,k,l)')/nbrOfRealizations;
              
        end
    end
end


Gp_MR = sum(Gp_MR_total,4);
%------------Calculation of LSFD coefficients
for k = 1:K
    
    b = signal_MR(:,:,k);
    A(:,:,k) = p(k)*((Gp_MR(:,:,k)) + scaling_MR(:,:,k))\b;
        
end


%Compute the SE
for k = 1:K
    
    b = signal_MR(:,:,k);
    numerator = sqrt(p(k))*A(:,:,k)'*b;
    denominator = A(:,:,k)'*(Gp_MR(:,:,k))*A(:,:,k) - numerator*numerator' + A(:,:,k)'*scaling_MR(:,:,k)*A(:,:,k);
    SE(k) = prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)));
        
end

        
