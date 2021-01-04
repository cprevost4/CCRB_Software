clear all
close all
clc
rng(107) %For reproducibility

%% Groundtruth data

%Dimensions oF CP model
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
A1 = P1*A2;
B1 = P2*B2;
C2 = Pm*C1;

%Groundtruth tensors
X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2}); 

%% Uncoupled CRB

for s=1:length(sigma_n2)
    [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2,B2,C2,sigma_n2(s));
    CRB_C2(s) = sum(diag(CRB2(1:dim2(3)*F,1:dim2(3)*F,s)));
    CRB_AB2(s) = sum(diag(CRB2(dim2(3)*F+1:end,dim2(3)*F+1:end,s)));
end

for s=1:length(sigma_n1)
    [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled2(A1,B1,C1,sigma_n1(s));
    CRB_C1(s) = sum(diag(CRB1(1:dim1(3)*F,1:dim1(3)*F,s)));
    CRB_AB1(s) = sum(diag(CRB1(dim1(3)*F+1:end,dim1(3)*F+1:end,s)));
    CRB_psi1(s) = CRB_AB2+CRB_C1(s);
    CRB_psi2(s) = CRB_AB1(s)+CRB_C2;
    FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2);
end

%% Uncoupled CRB for reconstructed SRI

[J_12,J_13,~] = vec_unfold([dim2(1)-1 dim2(2)-1 dim1(3)]);
S_A = kron(kr(C1,B2(2:end,:)),eye(dim2(1)-1));
S_B = J_12'*kron(kr(C1,A2(2:end,:)),eye(dim2(2)-1));
S_C = J_13'*kron(kr(B2(2:end,:),A2(2:end,:)),eye(dim1(3)));
H = [S_A S_B S_C];

for s=1:length(sigma_n1)
    [CRB(:,:,s)] = crb_allatonce(A1,B1,C1,A2,B2,C2,sigma_n1(s),sigma_n2);
    CRB_theta(:,:,s) = CRB(1:(dim1(3)+dim2(1)+dim2(2)-2)*F,1:(dim1(3)+dim2(1)+dim2(2)-2)*F,s);
    CRB_rec(:,:,s) = H*CRB_theta(:,:,s)*H';
    CRB_Y(s) = sum(diag(CRB_rec(:,:,s)));
end

%% CCRB - fully coupled model

H3 = [-kron(eye(F),Pm) zeros(dim2(3)*F,(dim1(1)+dim1(2))*F) eye(dim2(3)*F) zeros(dim2(3)*F,(dim2(1)+dim2(2)-2)*F); ...
        zeros((dim1(1)+dim1(2))*F,dim1(3)*F) eye((dim1(1)+dim1(2))*F) zeros((dim1(1)+dim1(2))*F,dim2(3)*F) -blkdiag(kron(eye(F),P1(:,2:end)),kron(eye(F),P2(:,2:end)))];
U = null(H3);

for s=1:length(sigma_n1)
    CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
    CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F,1:dim1(3)*F,s)));
    CCRB2_AB1(s) = sum(diag(CCRB(dim1(3)*F+1:(sum(dim1))*F,dim1(3)*F+1:(sum(dim1))*F,s)));
    CCRB2_C2(s) = sum(diag(CCRB((sum(dim1))*F+1:(dim2(3)+sum(dim1))*F,(sum(dim1))*F+1:(dim2(3)+sum(dim1))*F,s)));
    CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1))*F+1:end,(dim2(3)+sum(dim1))*F+1:end,s)));
    CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
    CCRB1_psi2(s) = CCRB2_AB1(s)+CCRB2_C2(s);
end

%% Reparametrized CRB

[J_12,J_13,~] = vec_unfold([dim2(1)-1 dim2(2)-1 dim1(3)]);
S_A = kron(kr(C1,B2(2:end,:)),eye(dim2(1)-1));
S_B = J_12'*kron(kr(C1,A2(2:end,:)),eye(dim2(2)-1));
S_C = J_13'*kron(kr(B2(2:end,:),A2(2:end,:)),eye(dim1(3)));
H = [S_A S_B S_C];
H2 = blkdiag(kron(eye(F),P1(:,2:end)),kron(eye(F),P2(:,2:end)),kron(eye(F),Pm));

for s=1:length(SNR1)
    CCRB_theta(:,:,s) = ccrb2(A2,B2,C1,dim1,dim2,sigma_n1(s),sigma_n2,P1,P2,Pm);
    CCRB_alpha(:,:,s) = H2*CCRB_theta(:,:,s)*H2';
    CCRB_rec(:,:,s) = H*CCRB_theta(:,:,s)*H'; 
    
    CCRB_AB2(s) = sum(diag(CCRB_theta(1:(dim2(1)+dim2(2)-2)*F,1:(dim2(1)+dim2(2)-2)*F,s)));
    CCRB_C1(s) = sum(diag(CCRB_theta((dim2(1)+dim2(2)-2)*F+1:end,(dim2(1)+dim2(2)-2)*F+1:end,s)));
    CCRB_AB1(s) = sum(diag(CCRB_alpha(1:(dim1(1)+dim1(2)-2)*F,1:(dim1(1)+dim1(2)-2)*F,s)));
    CCRB_C2(s) = sum(diag(CCRB_alpha((dim1(1)+dim1(2)-2)*F+1:end,(dim1(1)+dim1(2)-2)*F+1:end,s)));
    CCRB_Y(s) = sum(diag(CCRB_rec(:,:,s)));
end

%% Figures

load('performance2.mat')
%bounds_nl = load('bounds.mat');

figure
subplot(1,2,1)
semilogy(SNR1,CCRB1_psi1,'k--','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,mse_psi1_c,'bo','Linewidth',0.7,'MarkerSize',8); hold on
title('Performance bounds for $\widetilde{\mathbf{\theta}}$','Interpreter','latex')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
legend('$\textbf{CCRB}$','MSE - STEREO','Interpreter','latex')
xlim([5 60])
set(gca,'FontName','Times','FontSize',14)
subplot(1,2,2)
semilogy(SNR1,CCRB_Y,'k--','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,mse_Y_c,'bo','Linewidth',0.7,'MarkerSize',8); hold on
title('Performance bounds for $\mathbf{x}$','Interpreter','latex')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
legend('$\textbf{CCRB}$','MSE - STEREO','Interpreter','latex')
xlim([5 60])
set(gca,'FontName','Times','FontSize',14)
saveas(gcf,'figures/fig4.fig')

figure
semilogy(SNR1,CCRB1_psi1+CCRB1_psi2,'k--','Linewidth',1,'MarkerSize',8); hold on
semilogy(SNR1,CCRB_AB1+CCRB_AB2+CCRB_C1+CCRB_C2,'r+','Linewidth',1,'MarkerSize',8); hold on
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
legend('$\textbf{CCRB}$','Reparametrized CRB','Interpreter','latex')
xlim([5 60])
set(gca,'FontName','Times','FontSize',14)
saveas(gcf,'figures/fig1.fig')

