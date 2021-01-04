clear all
close all
clc
rng(107) %For reproducibility

%% Groundtruth data

%Dimensions of CP model
dim1 = [6,6,16]; dim2 = [18,18,8]; F = 3;

%Noise on tensors
SNR1 = 5:5:60; 
for s=1:length(SNR1)
    sigma_n1(s) = 10^(-SNR1(s)/10);
end
SNR2 = 20; sigma_n2 = 10^(-SNR2/10);
 
%Degradation matrices
q = 3; phi = gauss_kernel(q); phi = phi/norm(phi);
H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);
d = 3; S = eye(d*dim1(1)); S = S(1:d:end,:);
P1 = S*H; P2 = P1;
h = 2; phi2 = (1/h)*ones(1,h);
H2 = toeplitz([phi2(1),zeros(1,dim1(3)-1)], [phi2 zeros(1,dim1(3)-h)]);
d = 2; S2 = eye(d*dim2(3)); S2 = S2(1:d:end,:);
Pm = S2*H2;

%CP factors
A2 = randn(dim2(1),F); B2 = randn(dim2(2),F); C1 = randn(dim1(3),F);
A2(1,:) = 1; B2(1,:) = 1;
alpha = P1(1,:)*A2; beta = P2(1,:)*B2;
A1 = P1*A2*diag(1./alpha);
B1 = P2*B2*diag(1./beta);
C2 = Pm*C1*diag(1./(alpha.*beta));

%Groundtruth tensors
X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2}); 

%% Uncoupled CRB

for s=1:length(sigma_n2)
    [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2,B2,C2,sigma_n2(s));
    CRB_C2(s) = sum(diag(CRB2(1:dim2(3)*F,1:dim2(3)*F,s)));
    CRB_AB2(s) = sum(diag(CRB2(dim2(3)*F+1:end,dim2(3)*F+1:end,s)));
end

for s=1:length(sigma_n1)
    [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled(A1,B1,C1,sigma_n1(s));
    CRB_C1(s) = sum(diag(CRB1(1:dim1(3)*F,1:dim1(3)*F,s)));
    CRB_AB1(s) = sum(diag(CRB1(dim1(3)*F+1:end,dim1(3)*F+1:end,s)));
    CRB_psi1(s) = CRB_AB2+CRB_C1(s);
    CRB_psi2(s) = CRB_AB1(s)+CRB_C2;
    FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2); 
end

%% Uncoupled CRB for reconstructed tensor

%Transformation matrix
[J_12,J_13,~] = vec_unfold([dim2(1)-1 dim2(2)-1 dim1(3)]);
S_A = kron(kr(C1,B2(2:end,:)),eye(dim2(1)-1));
S_B = J_12'*kron(kr(C1,A2(2:end,:)),eye(dim2(2)-1));
S_C = J_13'*kron(kr(B2(2:end,:),A2(2:end,:)),eye(dim1(3)));
H = [S_A S_B S_C];

% Uncoupled CRB for Psi
for s=1:length(sigma_n1)
    [CRB(:,:,s)] = crb_allatonce(A1,B1,C1,A2,B2,C2,sigma_n1(s),sigma_n2);
    CRB_theta(:,:,s) = CRB(1:(dim1(3)+dim2(1)+dim2(2)-2)*F,1:(dim1(3)+dim2(1)+dim2(2)-2)*F,s);
    CRB_rec(:,:,s) = H*CRB_theta(:,:,s)*H';

    CRB_Y(s) = sum(diag(CRB_rec(:,:,s)));
end

%% CCRB - fully coupled model

%Jacobian of constraints
[~, H3] = transformation_jac(A2,B2,C1,P1,P2,Pm,alpha,beta,dim1,dim2);
U = null(H3);

%Standard CCRB
for s=1:length(sigma_n1)
    CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
    CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F,1:dim1(3)*F,s)));
    CCRB2_AB1(s) = sum(diag(CCRB(dim1(3)*F+1:(sum(dim1)-2)*F,dim1(3)*F+1:(sum(dim1)-2)*F,s)));
    CCRB2_C2(s) = sum(diag(CCRB((sum(dim1)-2)*F+1:(dim2(3)+sum(dim1)-2)*F,(sum(dim1)-2)*F+1:(dim2(3)+sum(dim1)-2)*F,s)));
    CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*F+1:end,(dim2(3)+sum(dim1)-2)*F+1:end,s)));
    CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
    CCRB1_psi2(s) = CCRB2_AB1(s)+CCRB2_C2(s);
end

%% Reparametrized CRB

%Transformation matrix
[H2, ~] = transformation_jac(A2,B2,C1,P1,P2,Pm,alpha,beta,dim1,dim2);
[J_12,J_13,~] = vec_unfold([dim2(1)-1 dim2(2)-1 dim1(3)]);
S_A = kron(kr(C1,B2(2:end,:)),eye(dim2(1)-1));
S_B = J_12'*kron(kr(C1,A2(2:end,:)),eye(dim2(2)-1));
S_C = J_13'*kron(kr(B2(2:end,:),A2(2:end,:)),eye(dim1(3)));
H = [S_A S_B S_C];

%Reparametrized CRB
for s=1:length(SNR1)
    CCRB_theta(:,:,s) = ccrb3(A2,B2,C1,dim1,dim2,sigma_n1(s),sigma_n2,P1,P2,Pm);
    CCRB_alpha(:,:,s) = H2*CCRB_theta(:,:,s)*H2';
    CCRB_rec(:,:,s) = H*CCRB_theta(:,:,s)*H'; 
    
    CCRB_AB2(s) = sum(diag(CCRB_theta(1:(dim2(1)+dim2(2)-2)*F,1:(dim2(1)+dim2(2)-2)*F,s)));
    CCRB_C1(s) = sum(diag(CCRB_theta((dim2(1)+dim2(2)-2)*F+1:end,(dim2(1)+dim2(2)-2)*F+1:end,s)));
    CCRB_AB1(s) = sum(diag(CCRB_alpha(1:(dim1(1)+dim1(2)-2)*F,1:(dim1(1)+dim1(2)-2)*F,s)));
    CCRB_C2(s) = sum(diag(CCRB_alpha((dim1(1)+dim1(2)-2)*F+1:end,(dim1(1)+dim1(2)-2)*F+1:end,s)));
    CCRB_Y(s) = sum(diag(CCRB_rec(:,:,s)));
end

%% Standard Blind-CCRB

%Jacobian
dC2_dC1 = kron(diag(1./(alpha.*beta)),Pm);
dC2_dA2 = -kron(diag(1./(beta.*(alpha.^2))),Pm)*kr(eye(F),C1)*kron(eye(F),P1(1,2:end));
dC2_dB2 = -kron(diag(1./(alpha.*(beta.^2))),Pm)*kr(eye(F),C1)*kron(eye(F),P2(1,2:end));
H3 = [-dC2_dC1 zeros(dim2(3)*F,(dim1(1)+dim1(2)-2)*F) eye(dim2(3)*F) -dC2_dA2 -dC2_dB2];
U = null(H3);

%Blind-CCRB
for s=1:length(sigma_n1)
    CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
    CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F,1:dim1(3)*F,s)));
    CCRB2_AB1(s) = sum(diag(CCRB(dim1(3)*F+1:(sum(dim1)-2)*F,dim1(3)*F+1:(sum(dim1)-2)*F,s)));
    CCRB2_C2(s) = sum(diag(CCRB((sum(dim1)-2)*F+1:(dim2(3)+sum(dim1)-2)*F,(sum(dim1)-2)*F+1:(dim2(3)+sum(dim1)-2)*F,s)));
    CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*F+1:end,(dim2(3)+sum(dim1)-2)*F+1:end,s)));
    CCRB2_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
    CCRB2_psi2(s) = CCRB2_AB1(s)+CCRB2_C2(s);
end

%% Reparametrized Blind-CCRB

H3 = [zeros(dim2(3)*F,(dim1(1)+dim1(2)-2)*F) dC2_dA2 dC2_dB2 dC2_dC1];

for s=1:length(SNR1)
    BlindCCRB(:,:,s) = ccrb4(A1,A2,B1,B2,C1,P1,P2,alpha,beta,dim1,dim2,sigma_n1(s),sigma_n2,Pm);
    BCCRB(:,:,s) = H3*BlindCCRB(:,:,s)*H3';
    
    CCRB3(s) = sum(diag(BlindCCRB(:,:,s))) + sum(diag(BCCRB(:,:,s)));
    CCRB3_psi1(s) = sum(diag(BlindCCRB(:,:,s)));
end

%% Blind-CCRB for reconstructed tensor

for s=1:length(sigma_n1)
    CCRB4(:,:,s) = [CCRB((sum(dim1)-2)*F+1:end,(sum(dim1)-2)*F+1:end,s) CCRB((sum(dim1)-2)*F+1:end,1:(sum(dim1)-2)*F,s);
                    CCRB(1:(sum(dim1)-2)*F,(sum(dim1)-2)*F+1:end,s) CCRB(1:(sum(dim1)-2)*F,1:(sum(dim1)-2)*F,s)];
    CCRB_rec_b(:,:,s) = H*CCRB4(dim2(3)*F+1:(sum(dim2)+dim1(3)-2)*F,dim2(3)*F+1:(sum(dim2)+dim1(3)-2)*F,s)*H';
    CCRB_Y_b(s) = sum(diag(CCRB_rec_b(:,:,s)));

    test(s) = sum(diag(CCRB4(dim2(3)*F+1:(sum(dim2)+dim1(3)-2)*F,dim2(3)*F+1:(sum(dim2)+dim1(3)-2)*F,s)));
end

%% Figures

load('performance.mat')

figure
subplot(1,2,1)
semilogy(SNR1,CRB_psi1,'k-','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,mse_psi1,'r+','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,CCRB1_psi1,'k--','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,mse_psi1_c,'bo','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,CCRB2_psi1,'k:','Linewidth',2,'MarkerSize',8); hold on
semilogy(SNR1,mse_psi1_b,'gd','Linewidth',1,'MarkerSize',8); hold on
title('Performance bounds for $\widetilde{\mathbf{\theta}}$','Interpreter','latex')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
legend('$\textbf{CRB}$','MSE - Uncoupled ALS','$\textbf{CCRB}$','MSE - STEREO','$\textbf{Blind-CCRB}$','MSE - Blind-STEREO','Interpreter','latex')
set(gca,'FontName','Times','FontSize',14)
xlim([5 60])
%ylim([2e-4 2])
subplot(1,2,2)
semilogy(SNR1,CCRB_Y,'k--','Linewidth',1,'MarkerSize',8); hold on
%semilogy(SNR1,mse_Y_u,'r+','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,CCRB_Y_b,'k:','Linewidth',2,'MarkerSize',8); hold on
semilogy(SNR1,mse_Y_b,'gd','Linewidth',1,'MarkerSize',8); hold on
%semilogy(SNR1,CRB_Y,'k-','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,mse_Y_c,'bo','Linewidth',1,'MarkerSize',8)
title('Performance bounds for $\mathbf{x}$','Interpreter','latex')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',14)
xlim([5 60])
saveas(gcf,'figures/fig5.fig')


figure
semilogy(SNR1,CRB_psi1+CRB_psi2,'k-','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,CCRB1_psi1+CCRB1_psi2,'k--','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,CCRB_AB2+CCRB_C1+CCRB_AB1+CCRB_C2,'r+','Linewidth',0.7,'MarkerSize',8); hold on
semilogy(SNR1,CCRB2_psi1+CCRB2_psi2,'k:','Linewidth',2,'MarkerSize',8); hold on
semilogy(SNR1,CCRB3,'bo','Linewidth',0.7,'MarkerSize',8); hold on
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
legend('$\textbf{CRB}$','$\textbf{CCRB}$','Reparametrized CRB','$\textbf{Blind-CCRB}$','Blind reparametrized CRB','Interpreter','latex')
xlim([5 60]); %ylim([5e-4,2])
set(gca,'FontName','Times','FontSize',14)
saveas(gcf,'figures/fig2.fig')


