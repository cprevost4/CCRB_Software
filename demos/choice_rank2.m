clear all
close all
clc

%% Generate model

% Dimensions
F_th = 3; dim1 = [6,6,16]; dim2 = [18,18,8];
rng(107)

% Degradation matrices
q = 3; phi = gauss_kernel(q);
H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);
d = 3; S = eye(d*dim1(1)); S = S(1:d:end,:);
P1 = S*H; P2 = P1;
Pm = eye(dim1(3));
Pm(1:end-1,2:end) = Pm(1:end-1,2:end) + eye(dim1(3)-1); Pm = Pm(1:2:end,:);
Pm = Pm/2;

% True factors
A2_th = randn(dim2(1),F_th); B2_th = randn(dim2(2),F_th); C1_th = randn(dim1(3),F_th);
A2_th(1,:) = 1; B2_th(1,:) = 1;
A1_th = P1*A2_th; B1_th = P2*B2_th; C2_th = Pm*C1_th;


% SNR
SNR1 = 5:5:60; 
for s=1:length(SNR1)
    sigma_n1(s) = 10^(-SNR1(s)/10);
end
SNR2 = 20; sigma_n2 = 10^(-SNR2/10); 

%% Run bounds

Nreal = 100;

for n=1:Nreal
    
    n
    
    % Additional columns
    A2_add = randn(dim2(1),13); B2_add = randn(dim2(2),13); C1_add = randn(dim1(3),13);
    A2_add(1,:) = 1; B2_add(1,:) = 1;
    A1_add = P1*A2_add; B1_add = P2*B2_add; C2_add = Pm*C1_add;

    for F=3:16

        clear CRB2 FIM2 CRB1 FIM1 FIM CCRB

        A1 = [A1_th A1_add(:,1:F-F_th)]; B1 = [B1_th B1_add(:,1:F-F_th)]; C1 = [C1_th C1_add(:,1:F-F_th)]; 
        A2 = [A2_th A2_add(:,1:F-F_th)]; B2 = [B2_th B2_add(:,1:F-F_th)]; C2 = [C2_th C2_add(:,1:F-F_th)]; 

        for s=1:length(sigma_n2)
            [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2,B2,C2,sigma_n2(s));
        end
        for s=1:length(sigma_n1)
            [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled2(A1,B1,C1,sigma_n1(s));
            FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2);
        end

        H3 = [-kron(eye(F),Pm) zeros(dim2(3)*F,(dim1(1)+dim1(2))*F) eye(dim2(3)*F) zeros(dim2(3)*F,(dim2(1)+dim2(2)-2)*F); ...
                zeros((dim1(1)+dim1(2))*F,dim1(3)*F) eye((dim1(1)+dim1(2))*F) zeros((dim1(1)+dim1(2))*F,dim2(3)*F) -blkdiag(kron(eye(F),P1(:,2:end)),kron(eye(F),P2(:,2:end)))];
        U = null(H3);

        for s=1:length(sigma_n1)
            CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';

            CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F_th,1:dim1(3)*F_th,s)));
            CCRB2_A2(s) = sum(diag(CCRB((sum(dim1))*F+1:(sum(dim1))*F+(dim2(1)-1)*F_th,(sum(dim1))*F+1:(sum(dim1))*F+(dim2(1)-1)*F_th,s)));
            CCRB2_B2(s) = sum(diag(CCRB((sum(dim1)+dim2(1)-1)*F+1:(sum(dim1)+dim2(1)-1)*F+(dim2(2)-1)*F_th,(sum(dim1)+dim2(1)-1)*F+1:(sum(dim1)+dim2(1)-1)*F+(dim2(2)-1)*F_th,s)));
            CCRB1_psi1(s) = CCRB2_A2(s)+CCRB2_B2(s)+CCRB2_C1(s);
        end

        bound(F-2,n,:) = CCRB1_psi1;
    end

end

%% Figure

figure(2)
surf(5:5:60,3:16,bound); hold on
set(gca,'zscale','log')
set(gca,'ColorScale','log')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
ylabel('Tensor rank $N$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',14)
xlim([5 60]); ylim([3 16])
colorbar
saveas(gcf,'figures/fig3.fig')