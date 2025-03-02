%%=============================================================
%The file is a test of simulation codes of the paper:
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


warning('off');

clear all
close all


M = 20;
K = 10;
L = 2;
N = 2;
tau_c = 200;
nbrOfSetups = 2; %Number of realizations for AP/UE locations
nbrOfRealizations = 800; %Number of channel realizations for each AP/UE location
w = 2; %Pilot Reuse factor
p = 0.2; %200 mW
%Create the power vector for all UEs (The uplink power is the same
%(p)at each UE)
p_u = p*ones(1,K);


SE_Monte_MMSE_Combining_Level2 = zeros(K,nbrOfSetups,length(N));
SE_Theoretical_MR_Combining_Level2 = zeros(K,nbrOfSetups,length(N));
SE_Monte_MR_Combining_Level2 = zeros(K,nbrOfSetups,length(N));
% 
SE_Monte_MMSE_Combining_Level3 = zeros(K,nbrOfSetups,length(N));
SE_Monte_MR_Combining_Level3 = zeros(K,nbrOfSetups,length(N));
SE_Theoretical_MR_Combining_Level3 = zeros(K,nbrOfSetups,length(N));

SE_Monte_MMSE_Combining_Level1 = zeros(K,nbrOfSetups,length(N));
SE_Monte_MMSE_Combining_Level4 = zeros(K,nbrOfSetups,length(N));

SE_Monte_MR_Combining_Level1 = zeros(K,nbrOfSetups,length(N));
SE_Monte_MR_Combining_Level4 = zeros(K,nbrOfSetups,length(N));



for n = 1:length(N)

    tau_p = K*N(n)/w;
    p_antenna = p_u/N(n);

    for i = 1:nbrOfSetups
        
        [channelGain] = RandomAP_generateSetup_Rician_Multi_Antenna(M,K,1,1);
        [HH,R_Vec,H_Vec,Omega_couple] = functionChannelGeneration(channelGain,M,K,N(n),L,nbrOfRealizations);
        disp(['ChannelGeneration of setup ' num2str(i)]);
        
        [Pset] = functionPilotAllocation( R_Vec,M,K,L*N(n),tau_p/N(n),p_antenna);
        % disp(['PilotAllocation of setup ' num2str(i)]);
        
        [Hhat_MMSE] = functionChannelEstimates_MMSE(H_Vec,R_Vec,nbrOfRealizations,M,K,L,N(n),tau_p,p_antenna,Pset);
        
        [Phi,Omega,C_MMSE,C_MMSE_MMSE_Combining] = functionMatrixGeneration(R_Vec,M,K,L,N(n),tau_p,p_antenna,Pset);
        disp(['MatrixGeneration of setup ' num2str(i)]);
        
        [V_MMSE_Combining] = functionCompute_MMSE_Combining_Matrix(Hhat_MMSE,C_MMSE_MMSE_Combining,nbrOfRealizations,L,N(n),K,M,p_antenna);
        disp(['L-MMSE combining matrix of setup ' num2str(i)]);
        
        %--------Level3
        %----Monte
        [SE_MR_Combining_LSFD] = functionComputeMonteCarloSE_UL_LSFD(Hhat_MMSE,HH,tau_c,tau_p,nbrOfRealizations,N(n),L,K,M,p_antenna);
        disp(['MR Combining Monte Compution of Level3 ' num2str(i)]);
        [SE_MMSE_Combining_LSFD] = functionComputeMonteCarloSE_UL_LSFD(V_MMSE_Combining,HH,tau_c,tau_p,nbrOfRealizations,N(n),L,K,M,p_antenna);
        disp(['MMSE Combining Monte Compution of Level3 ' num2str(i)]);
       
        %----Theoretical
        [SE_MR_th_LSFD] = functionTheoreticalCellFreeULSE_MMSE_LSFD( R_Vec,Phi,Omega,M,K,L,N(n),p_antenna,tau_p,tau_c,Pset); 
        disp(['MR Combining Theoretical SE of Level3 ' num2str(i)]);
        
        %--------Level2
        %----Monte       
        [SE_MR_Combining_Level2]= functionComputeMonteCarloSE_UL_L2(Hhat_MMSE,HH,tau_c,tau_p,nbrOfRealizations,M,K,L,N(n),p_antenna);
        disp(['MR Combining Monte Compution of Level2 ' num2str(i)]);
        [SE_MMSE_Combining_Level2]= functionComputeMonteCarloSE_UL_L2(V_MMSE_Combining,HH,tau_c,tau_p,nbrOfRealizations,M,K,L,N(n),p_antenna);
        disp(['MMSE Combining Monte Compution of Level2 ' num2str(i)]);
        
        %----Theoretical
        [SE_MR_th_Level2] = functionTheoreticalCellFreeULSE_MMSE_Level2(R_Vec,Phi,Omega,M,K,L,N(n),p_antenna,tau_p,tau_c,Pset);
        disp(['MR Combining Theoretical SE of Level2 ' num2str(i)]);
        
        %--------Level4 & Level1
        [SE_MR_Combining_Level4,SE_MR_Combining_Level1] = functionComputeSE_Fully_Centralized_Small_Cell(Hhat_MMSE,Hhat_MMSE,C_MMSE_MMSE_Combining,tau_c,tau_p,nbrOfRealizations,N(n),L,K,M,p_antenna,1);
        disp(['SE with MR combining vector of L1 and L4 ' num2str(i)]);

        [SE_MMSE_Combining_Level4,SE_MMSE_Combining_Level1] = functionComputeSE_Fully_Centralized_Small_Cell(V_MMSE_Combining,Hhat_MMSE,C_MMSE_MMSE_Combining,tau_c,tau_p,nbrOfRealizations,N(n),L,K,M,p_antenna,0);
        disp(['SE with MMSE combining vector of L1 and L4 ' num2str(i)]);
   
        SE_Monte_MMSE_Combining_Level2(:,i,n) = SE_MMSE_Combining_Level2;
        SE_Theoretical_MR_Combining_Level2(:,i,n) = SE_MR_th_Level2;
        SE_Monte_MR_Combining_Level2(:,i,n) = SE_MR_Combining_Level2;

        SE_Monte_MMSE_Combining_Level3(:,i,n) = SE_MMSE_Combining_LSFD;
        SE_Monte_MR_Combining_Level3(:,i,n) = SE_MR_Combining_LSFD;
        SE_Theoretical_MR_Combining_Level3(:,i,n) = SE_MR_th_LSFD;   
        
        SE_Monte_MMSE_Combining_Level1(:,i,n) = SE_MMSE_Combining_Level1;
        SE_Monte_MMSE_Combining_Level4(:,i,n) = SE_MMSE_Combining_Level4;

        SE_Monte_MR_Combining_Level1(:,i,n) = SE_MR_Combining_Level1;
        SE_Monte_MR_Combining_Level4(:,i,n) = SE_MR_Combining_Level4;

        
    end
    disp([num2str(n) ' UE-Antenna out of ' num2str(N)]);
end

