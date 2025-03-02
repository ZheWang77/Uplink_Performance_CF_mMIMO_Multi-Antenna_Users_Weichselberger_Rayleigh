function [SE_CC] = functionTheoreticalCellFreeULSE_MMSE_Level2( R_Vec,Phi,Omega,M,K,L,N,p,tau_p,tau_c,Pset)
%%=============================================================
%This function is used to compute the theoretical uplink SE for Level 2
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
Zk = zeros(N,N,K,M);     %Zk_LMMSE
X_p1 = zeros(N,N,K,K,M); %Term Gama(1)
X_p2 = zeros(N,N,K,K,M); %Term Gama(2)
X_p3_1 = zeros(N,N,K,K,M); %Term m unequal m'
X_p3_2 = zeros(N,N,K,K,M);%Term m unequal m'
CCterm1 = zeros(N,N,K);  %D_k
CCterm2_P1_MMSE =  zeros(N,N,K,K);
CCterm3 =  zeros(N,N,K);

SE_CC = zeros(K,1);%Store the result

for k = 1:K
    for m = 1:M
        for n = 1:N
            for nn = 1:N
                
                Zk(n,nn,k,m) = trace(Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k));
                
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
%---------------------------------------------------------------------------------
X_p1_total = sum(X_p1,5);
X_p2_total = sum(X_p2,5);
X_p3_1_total = sum(X_p3_1,5);
X_p3_2_total = sum(X_p3_2,5);
Zk_total = sum(Zk,4);

for k = 1:K
    for l = 1:K
        
        CCterm1(:,:,k) = sqrt(p(k))*Zk_total(:,:,k);
        CCterm3(:,:,k) = Zk_total(:,:,k);
        CCterm2_P1_MMSE(:,:,k,l) = p(l)*X_p1_total(:,:,k,l);
        
        if any(l == Pset(:,k))
            
            X_p3_mm = zeros(N,N);
            
            for mm = 1:M
                
                X_p3_mm = X_p3_mm + X_p3_1(:,:,k,l,mm)*X_p3_2(:,:,k,l,mm);
                
            end
            
            CCterm2_P1_MMSE(:,:,k,l) = CCterm2_P1_MMSE(:,:,k,l) + p(l)*X_p2_total(:,:,k,l) - p(l)*X_p1_total(:,:,k,l) + p(l)*X_p3_1_total(:,:,k,l)*X_p3_2_total(:,:,k,l) - p(l)*X_p3_mm;
            clear X_p3_mm
   
        end
    end
end

CCterm2 = sum(CCterm2_P1_MMSE,4);


for k = 1:K
    
    SE_CC(k) = prelogFactor*real(log2(det(eye(N) + CCterm1(:,:,k)'/((CCterm2(:,:,k)) - CCterm1(:,:,k)*CCterm1(:,:,k)' + CCterm3(:,:,k))*CCterm1(:,:,k))));
    
end
    

