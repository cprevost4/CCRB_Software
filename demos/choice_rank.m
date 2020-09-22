clear all
close all
clc
rng(107) %For reproducibility

%% For scenario 2

%Generate data
F_th = 3; dim1 = [6,6,16]; dim2 = [18,18,8];

q = 3; phi = gauss_kernel(q);
H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);
d = 3; S = eye(d*dim1(1)); S = S(1:d:end,:);
P1 = S*H; P2 = P1;
Pm = eye(dim1(3));
Pm(1:end-1,2:end) = Pm(1:end-1,2:end) + eye(dim1(3)-1); Pm = Pm(1:2:end,:);
Pm = Pm/2;

A2_th = randn(dim2(1),F_th); B2_th = randn(dim2(2),F_th); C1_th = randn(dim1(3),F_th);
A2_th(1,:) = 1; B2_th(1,:) = 1;

alpha = P1(1,:)*A2_th; beta = P2(1,:)*B2_th;
A1_th = P1*A2_th*diag(1./alpha);
B1_th = P2*B2_th*diag(1./beta);
C2_th = Pm*C1_th*diag(1./(alpha.*beta));
A1_add = randn(dim1(1),15); B1_add = randn(dim1(2),15); C1_add = randn(dim1(3),15);
A2_add = randn(dim2(1),15); B2_add = randn(dim2(2),15); C2_add = randn(dim2(3),15);
A1_add(1,:) = 1; B1_add(1,:) = 1; A2_add(1,:) = 1; B2_add(1,:) = 1;

SNR1 = 5:5:60; 
for s=1:length(SNR1)
    sigma_n1(s) = F_th*10^(-SNR1(s)/10);
end
SNR2 = 20; sigma_n2 = F_th*10^(-SNR2/10);

%Compute bounds
for s=1:length(sigma_n2)
    [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2_th,B2_th,C2_th,sigma_n2(s));
end
for s=1:length(sigma_n1)
    [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled(A1_th,B1_th,C1_th,sigma_n1(s));
    FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2);
end

[~, H3] = transformation_jac(A2_th,B2_th,C1_th,P1,P2,Pm,alpha,beta,dim1,dim2);
U = null(H3);
for s=1:length(sigma_n1)
    CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
    CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F_th,1:dim1(3)*F_th,s)));
    CCRB2_AB1(s) = sum(diag(CCRB(dim1(3)*F_th+1:(sum(dim1)-2)*F_th,dim1(3)*F_th+1:(sum(dim1)-2)*F_th,s)));
    CCRB2_C2(s) = sum(diag(CCRB((sum(dim1)-2)*F_th+1:(dim2(3)+sum(dim1)-2)*F_th,(sum(dim1)-2)*F_th+1:(dim2(3)+sum(dim1)-2)*F_th,s)));
    CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*F_th+1:end,(dim2(3)+sum(dim1)-2)*F_th+1:end,s)));
    CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
    CCRB1_psi2(s) = CCRB2_AB1(s)+CCRB2_C2(s);
end
bound(1,:) = CCRB1_psi1;

for F=4:16
   
    SNR1 = 5:5:60; %Noise on first tensor
    
    for s=1:length(SNR1)
        sigma_n1(s) = F*10^(-SNR1(s)/10);
    end
    SNR2 = 20; sigma_n2 = F*10^(-SNR2/10);

    A1 = [A1_th A1_add(:,1:F-F_th)]; B1 = [B1_th B1_add(:,1:F-F_th)]; C1 = [C1_th C1_add(:,1:F-F_th)]; 
    A2 = [A2_th A2_add(:,1:F-F_th)]; B2 = [B2_th B2_add(:,1:F-F_th)]; C2 = [C2_th C2_add(:,1:F-F_th)]; 
    X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2});
    Y = cpdgen({A2(2:end,:),B2(2:end,:),C1});

    for s=1:length(sigma_n2)
        [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2(:,1:F_th),B2(:,1:F_th),C2(:,1:F_th),sigma_n2(s));
    end
    for s=1:length(sigma_n1)
        [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled(A1(:,1:F_th),B1(:,1:F_th),C1(:,1:F_th),sigma_n1(s));
        FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2);
    end

    [~, H3] = transformation_jac(A2(:,1:F_th),B2(:,1:F_th),C1(:,1:F_th),P1,P2,Pm,alpha(1:F_th),beta(1:F_th),dim1,dim2);
    U = null(H3);

    for s=1:length(sigma_n1)
        CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';

        CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F_th,1:dim1(3)*F_th,s)));
        CCRB2_AB1(s) = sum(diag(CCRB(dim1(3)*F_th+1:(sum(dim1)-2)*F_th,dim1(3)*F_th+1:(sum(dim1)-2)*F_th,s)));
        CCRB2_C2(s) = sum(diag(CCRB((sum(dim1)-2)*F_th+1:(dim2(3)+sum(dim1)-2)*F_th,(sum(dim1)-2)*F_th+1:(dim2(3)+sum(dim1)-2)*F_th,s)));
        CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*F_th+1:end,(dim2(3)+sum(dim1)-2)*F_th+1:end,s)));

        CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
        CCRB1_psi2(s) = CCRB2_AB1(s)+CCRB2_C2(s);
    end
    
    bound(F-2,:) = CCRB1_psi1;
end

figure(1)
surf(5:5:60,3:16,bound); hold on
set(gca,'zscale','log')
set(gca,'ColorScale','log')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
ylabel('Tensor rank $N$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',14)
colorbar
xlim([5 60]); ylim([2 16])
saveas(gcf,'figures/fig7.fig')

%% For scenario 1

clear all

F_th = 3; dim1 = [6,6,16]; dim2 = [18,18,8];

q = 3; phi = gauss_kernel(q);
H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);
d = 3; S = eye(d*dim1(1)); S = S(1:d:end,:);
P1 = S*H; P2 = P1;
Pm = eye(dim1(3));
Pm(1:end-1,2:end) = Pm(1:end-1,2:end) + eye(dim1(3)-1); Pm = Pm(1:2:end,:);
Pm = Pm/2;


A2_th = randn(dim2(1),F_th); B2_th = randn(dim2(2),F_th); C1_th = randn(dim1(3),F_th);
A2_th(1,:) = 1; B2_th(1,:) = 1;
alpha = P1(1,:)*A2_th; beta = P2(1,:)*B2_th;
A1_th = P1*A2_th;
B1_th = P2*B2_th;
C2_th = Pm*C1_th;
A1_add = randn(dim1(1),13); B1_add = randn(dim1(2),13); C1_add = randn(dim1(3),13);
A2_add = randn(dim2(1),13); B2_add = randn(dim2(2),13); C2_add = randn(dim2(3),13);
A2_add(1,:) = 1; B2_add(1,:) = 1;

SNR1 = 5:5:60; 
for s=1:length(SNR1)
    sigma_n1(s) = F_th*10^(-SNR1(s)/10);
end
SNR2 = 20; sigma_n2 = F_th*10^(-SNR2/10); 

for s=1:length(sigma_n2)
    [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2_th,B2_th,C2_th,sigma_n2(s));
end
for s=1:length(sigma_n1)
    [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled2(A1_th,B1_th,C1_th,sigma_n1(s));
    FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2); 
end

H3 = [-kron(eye(3),Pm) zeros(dim2(3)*3,(dim1(1)+dim1(2))*3) eye(dim2(3)*3) zeros(dim2(3)*3,(dim2(1)+dim2(2)-2)*3); ...
        zeros((dim1(1)+dim1(2))*3,dim1(3)*3) eye((dim1(1)+dim1(2))*3) zeros((dim1(1)+dim1(2))*3,dim2(3)*3) -blkdiag(kron(eye(3),P1(:,2:end)),kron(eye(3),P2(:,2:end)))];
U = null(H3); 
for s=1:length(sigma_n1)
    CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
    CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*3,1:dim1(3)*3,s)));
    CCRB2_AB1(s) = sum(diag(CCRB(dim1(3)*3+1:(sum(dim1))*3,dim1(3)*3+1:(sum(dim1))*3,s)));
    CCRB2_C2(s) = sum(diag(CCRB((sum(dim1))*3+1:(dim2(3)+sum(dim1))*3,(sum(dim1))*3+1:(dim2(3)+sum(dim1))*3,s)));
    CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*3+1:end,(dim2(3)+sum(dim1)-2)*3+1:end,s)));
    CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
    CCRB1_psi2(s) = CCRB2_AB1(s)+CCRB2_C2(s);
end

bound(1,:) = CCRB1_psi1;

for F=4:16
   
    SNR1 = 5:5:60;
    for s=1:length(SNR1)
        sigma_n1(s) = F*10^(-SNR1(s)/10);
    end
    SNR2 = 20; sigma_n2 = F*10^(-SNR2/10);

    A1 = [A1_th A1_add(:,1:F-F_th)]; B1 = [B1_th B1_add(:,1:F-F_th)]; C1 = [C1_th C1_add(:,1:F-F_th)]; 
    A2 = [A2_th A2_add(:,1:F-F_th)]; B2 = [B2_th B2_add(:,1:F-F_th)]; C2 = [C2_th C2_add(:,1:F-F_th)]; 

    X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2});
    Y = cpdgen({A2(2:end,:),B2(2:end,:),C1});

    for s=1:length(sigma_n2)
        [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2(:,1:3),B2(:,1:3),C2(:,1:3),sigma_n2(s));
    end
    for s=1:length(sigma_n1)
        [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled2(A1(:,1:3),B1(:,1:3),C1(:,1:3),sigma_n1(s));
        FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2);
    end

    H3 = [-kron(eye(3),Pm) zeros(dim2(3)*3,(dim1(1)+dim1(2))*3) eye(dim2(3)*3) zeros(dim2(3)*3,(dim2(1)+dim2(2)-2)*3); ...
            zeros((dim1(1)+dim1(2))*3,dim1(3)*3) eye((dim1(1)+dim1(2))*3) zeros((dim1(1)+dim1(2))*3,dim2(3)*3) -blkdiag(kron(eye(3),P1(:,2:end)),kron(eye(3),P2(:,2:end)))];
    U = null(H3);

    for s=1:length(sigma_n1)
        CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';

        CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*3,1:dim1(3)*3,s)));
        CCRB2_AB1(s) = sum(diag(CCRB(dim1(3)*3+1:(sum(dim1))*3,dim1(3)*3+1:(sum(dim1))*3,s)));
        CCRB2_C2(s) = sum(diag(CCRB((sum(dim1))*3+1:(dim2(3)+sum(dim1))*3,(sum(dim1))*3+1:(dim2(3)+sum(dim1))*3,s)));
        CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*3+1:end,(dim2(3)+sum(dim1)-2)*3+1:end,s)));

        CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
        CCRB1_psi2(s) = CCRB2_AB1(s)+CCRB2_C2(s);
    end

    bound(F-2,:) = CCRB1_psi1;
end


figure(2)
surf(5:5:60,3:16,bound); hold on
set(gca,'zscale','log')
set(gca,'ColorScale','log')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
ylabel('Tensor rank $N$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',14)
xlim([5 60]); ylim([3 16])
colorbar
saveas(gcf,'figures/fig6.fig')

