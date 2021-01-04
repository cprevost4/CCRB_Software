clear all
close all
clc
rng(107) %For reproducibility

%% Model with non-linear constraints (scenario 2)

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

%CRB
for s=1:length(sigma_n2)
    [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2,B2,C2,sigma_n2(s));
end
for s=1:length(sigma_n1)
    [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled(A1,B1,C1,sigma_n1(s));
    FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2); 
end

%CCRB (non-linear constraints)
[~, H3] = transformation_jac(A2,B2,C1,P1,P2,Pm,alpha,beta,dim1,dim2);
U = null(H3);
for s=1:length(sigma_n1)
    CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
    CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F,1:dim1(3)*F,s)));
    CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*F+1:end,(dim2(3)+sum(dim1)-2)*F+1:end,s)));
    CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
end

%CCRB (linear constraints)
H4 = [-kron(eye(F),Pm) zeros(dim2(3)*F,(dim1(1)+dim1(2)-2)*F) eye(dim2(3)*F) zeros(dim2(3)*F,(dim2(1)+dim2(2)-2)*F); ...
        zeros((dim1(1)+dim1(2)-2)*F,dim1(3)*F) eye((dim1(1)+dim1(2)-2)*F) zeros((dim1(1)+dim1(2)-2)*F,dim2(3)*F) -blkdiag(kron(eye(F),P1(2:end,2:end)),kron(eye(F),P2(2:end,2:end)))];
U2 = null(H4);
for s=1:length(sigma_n1)
    CCRB2(:,:,s) = U2*inv(U2'*FIM(:,:,s)*U2)*U2';
    CCRB3_C1(s) = sum(diag(CCRB2(1:dim1(3)*F,1:dim1(3)*F,s)));
    CCRB3_AB2(s) = sum(diag(CCRB2((dim2(3)+sum(dim1)-2)*F+1:end,(dim2(3)+sum(dim1)-2)*F+1:end,s)));
    CCRB2_psi1(s) = CCRB3_AB2(s)+CCRB3_C1(s);
end

%Performance of STEREO with model (32)
load('performance.mat')

%Plot results
figure
semilogy(SNR1,CCRB1_psi1,'k--'); hold on; semilogy(SNR1,CCRB2_psi1,'r--+');
hold on; semilogy(SNR1,mse_psi1_c,'bo');
title('Fig. 6 (left)')
saveas(gcf,'figures/fig6_left.fig')


%% Model with linear constraints (scenario 1)

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
A1 = P1*A2;
B1 = P2*B2;
C2 = Pm*C1;
%Groundtruth tensors
X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2}); 

%CRB
for s=1:length(sigma_n2)
    [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2,B2,C2,sigma_n2(s));
end
for s=1:length(sigma_n1)
    [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled2(A1,B1,C1,sigma_n1(s));
    FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2);
end

%CCRB (linear constraints)
H3 = [-kron(eye(F),Pm) zeros(dim2(3)*F,(dim1(1)+dim1(2))*F) eye(dim2(3)*F) zeros(dim2(3)*F,(dim2(1)+dim2(2)-2)*F); ...
        zeros((dim1(1)+dim1(2))*F,dim1(3)*F) eye((dim1(1)+dim1(2))*F) zeros((dim1(1)+dim1(2))*F,dim2(3)*F) -blkdiag(kron(eye(F),P1(:,2:end)),kron(eye(F),P2(:,2:end)))];
U = null(H3);
for s=1:length(sigma_n1)
    CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
    CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F,1:dim1(3)*F,s)));
    CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1))*F+1:end,(dim2(3)+sum(dim1))*F+1:end,s)));
    CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
end

%CCRB (non-linear constraints)
alpha = P1(1,:)*A2; beta = P2(1,:)*B2;
dA1_dA2 = kron(diag(1./alpha),P1(:,2:end)) - ...
    (kron(diag(1./(alpha.^2)),P1(:,2:end))*kr(eye(F),A2(2:end,:))*kron(eye(F),P1(1,2:end)));
dB1_dB2 = kron(diag(1./beta),P2(:,2:end)) - ...
    (kron(diag(1./(beta.^2)),P2(:,2:end))*kr(eye(F),B2(2:end,:))*kron(eye(F),P2(1,2:end)));
dC2_dC1 = kron(diag(1./(alpha.*beta)),Pm);
dC2_dA2 = -kron(diag(1./(beta.*(alpha.^2))),Pm)*kr(eye(F),C1)*kron(eye(F),P1(1,2:end));
dC2_dB2 = -kron(diag(1./(alpha.*(beta.^2))),Pm)*kr(eye(F),C1)*kron(eye(F),P2(1,2:end));
H4 = [-dC2_dC1 zeros(dim2(3)*F,(dim1(1)+dim1(2))*F) eye(dim2(3)*F) -dC2_dA2 -dC2_dB2;
          zeros((dim1(1)+dim1(2))*F,dim1(3)*F) eye((dim1(1)+dim1(2))*F)...
                            zeros((dim1(1)+dim1(2))*F,dim2(3)*F) -blkdiag(dA1_dA2,dB1_dB2)];
U2 = null(H4);
for s=1:length(sigma_n1)
    CCRB2(:,:,s) = U2*inv(U2'*FIM(:,:,s)*U2)*U2';
    CCRB3_C1(s) = sum(diag(CCRB2(1:dim1(3)*F,1:dim1(3)*F,s)));
    CCRB3_AB2(s) = sum(diag(CCRB2((dim2(3)+sum(dim1))*F+1:end,(dim2(3)+sum(dim1))*F+1:end,s)));
    CCRB2_psi1(s) = CCRB3_AB2(s)+CCRB3_C1(s);
end

%Performance of STEREO with model (6)-(7)
load('performance2.mat')

%Plot results
figure
semilogy(SNR1,CCRB1_psi1,'r+--'); hold on; semilogy(SNR1,CCRB2_psi1,'k--');
hold on; semilogy(SNR1,mse_psi1_c,'bo');
title('Fig. 6 (right)')
saveas(gcf,'figures/fig6_right.fig')