function [SE] = functionTheoreticalCellFreeULSE_MMSE_LSFD( R_Vec,Phi,Omega,M,K,L,N,p,tau_p,tau_c,Pset)       
%%=============================================================
%This function is used to compute the theoretical uplink SE for the LSFD
%scheme with MR combining of the paper:
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

%Prepare to store the results
Gkk = zeros(M*N,N,K);     %Gkk
Gama_kl = zeros(M*N,M*N,K,K);
A = zeros(M*N,N,K);
S_k = zeros(M*N,M*N,K);
SE = zeros(K,1);

X_p1 = zeros(N,N,K,K,M); %Term Gama(1)
X_p2 = zeros(N,N,K,K,M); %Term Gama(2)
X_p3_1 = zeros(N,N,K,K,M); %Term m unequal m'
X_p3_2 = zeros(N,N,K,K,M);%Term m unequal m'

for k = 1:K
    for m = 1:M
        for n = 1:N
            for nn = 1:N
                
                Gkk((m-1)*N+n,nn,k) = trace(Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k));
                
                for l = 1:K  
                    for i = 1:N
                        
                        X_p1(n,nn,k,l,m) = X_p1(n,nn,k,l,m) + trace(R_Vec((i-1)*L+1:i*L,(i-1)*L+1:i*L,m,l)*Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k)); %Term Gama(1)
                        
                    end 
                    
                        if any(l == Pset(:,k))
                            
                            S_mk = R_Vec(:,:,m,k)/Phi(:,:,m,k);
                            F_mkl_1 = tau_p*S_mk*(Phi(:,:,m,k) - p(l)*tau_p*R_Vec(:,:,m,l))*S_mk';
                            F_mkl_2 = S_mk*R_Vec(:,:,m,l)*S_mk';
                            F_mkl_2_sqrt = F_mkl_2^(1/2);
                            R_ml_sqrt = R_Vec(:,:,m,l)^(1/2);

                            
                            if any(nn == n)
                                for i = 1:N
                                    for p1 = 1:N
                                        for p2 = 1:N
                                            
                                            X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) + p(k)*p(l)*tau_p^2*trace(F_mkl_2_sqrt((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*R_ml_sqrt((i-1)*L+1:i*L,(p1-1)*L+1:p1*L))*trace(F_mkl_2_sqrt((n-1)*L+1:n*L,(p2-1)*L+1:p2*L)*R_ml_sqrt((p2-1)*L+1:p2*L,(i-1)*L+1:i*L))...
                                                + p(k)*p(l)*tau_p^2*norm(F_mkl_2_sqrt((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*R_ml_sqrt((i-1)*L+1:i*L,(p2-1)*L+1:p2*L),'fro')^2;
                                        end
                                    end
                                    X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) + p(k)*trace(R_Vec((i-1)*L+1:i*L,(i-1)*L+1:i*L,m,l)*F_mkl_1((n-1)*L+1:n*L,(n-1)*L+1:n*L));
                                end
                                
                            else 
                                for i = 1:N
                                    for p1 = 1:N
                                        for p2 = 1:N
                                            
                                            X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) + p(k)*p(l)*tau_p^2*trace(F_mkl_2_sqrt((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*F_mkl_2_sqrt((nn-1)*L+1:nn*L,(p1-1)*L+1:p1*L)*R_ml_sqrt((i-1)*L+1:i*L,(p2-1)*L+1:p2*L)*R_ml_sqrt((p2-1)*L+1:p2*L,(i-1)*L+1:i*L))...
                                                + p(k)*p(l)*tau_p^2*trace(F_mkl_2_sqrt((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*R_ml_sqrt((i-1)*L+1:i*L,(p1-1)*L+1:p1*L))*trace(F_mkl_2_sqrt((nn-1)*L+1:nn*L,(p2-1)*L+1:p2*L)*R_ml_sqrt((p2-1)*L+1:p2*L,(i-1)*L+1:i*L));

                                        end
                                    end
                                    X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) + p(k)*trace(R_Vec((i-1)*L+1:i*L,(i-1)*L+1:i*L,m,l)*F_mkl_1((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                                end
                            end
                            
                            RR1 =  R_Vec(:,:,m,k)/Phi(:,:,m,k)*R_Vec(:,:,m,l);
                            RR2 =  R_Vec(:,:,m,l)/Phi(:,:,m,k)*R_Vec(:,:,m,k);
                            
                            X_p3_1(n,nn,k,l,m) = sqrt(p(l)*p(k))*tau_p*trace(RR1((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                            X_p3_2(n,nn,k,l,m) = sqrt(p(l)*p(k))*tau_p*trace(RR2((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                            
                            clear RR1 RR2 S_mk F_mkl_1 F_mkl_2 F_mkl_2_sqrt R_ml_sqrt
    
                        end
                end
            end
        end
    end
end

for k = 1:K
    for m = 1:M
        for mm = 1:M
            for l = 1:K
                if any(mm == m)
                    
                    Gama_kl((m-1)*N+1:m*N,(m-1)*N+1:m*N,k,l) = p(l)*X_p1(:,:,k,l,m);
                    
                end
                
                if any(l == Pset(:,k))
                    if any(mm == m)
                        
                        Gama_kl((m-1)*N+1:m*N,(m-1)*N+1:m*N,k,l) = p(l)*X_p2(:,:,k,l,m);
                    else
                        Gama_kl((m-1)*N+1:m*N,(mm-1)*N+1:mm*N,k,l) = p(l)*X_p3_1(:,:,k,l,m)*X_p3_2(:,:,k,l,mm);
                    end
                end
            end
        end
        S_k((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = Gkk((m-1)*N+1:m*N,:,k);
    end
end

                   
%------------Calculation of LSFD coefficients
Gama_kl_total = sum(Gama_kl,4);
for k = 1:K

    A(:,:,k) = p(k)*((Gama_kl_total(:,:,k)) + S_k(:,:,k))\Gkk(:,:,k);
    
end


%-----------Compute the SE
for k = 1:K
    
    numerator = sqrt(p(k))*A(:,:,k)'*Gkk(:,:,k);
    denominator = A(:,:,k)'*(Gama_kl_total(:,:,k))*A(:,:,k) - numerator*numerator' + A(:,:,k)'*S_k(:,:,k)*A(:,:,k);
    SE(k) = prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)));
    
end
                
