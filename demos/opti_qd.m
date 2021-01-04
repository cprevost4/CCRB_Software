clear all
close all
clc

%% Here you can run the full code to get the bounds, but it's long !
% Instead, the /data folder provides the bounds for ranks 3 and 16.

dim1 = [6,6,16]; dim2 = [36,36,8]; F = 16; %Dimensions of CP model
SNR1 = 5:5:60;
for s=1:length(SNR1)
    sigma_n1(s) = 10^(-SNR1(s)/10);
end
SNR2 = 20; sigma_n2 = 10^(-SNR2/10); %Noise on second tensor
Pm = eye(dim1(3));
Pm(1:end-1,2:end) = Pm(1:end-1,2:end) + eye(dim1(3)-1); Pm = Pm(1:2:end,:);
Pm = Pm/2;
rng(107)
A2 = randn(dim2(1),F); B2 = randn(dim2(2),F); C1 = randn(dim1(3),F);
A2(1,:) = 1; B2(1,:) = 1;
 
[J_12,J_13,~] = vec_unfold([dim2(1) dim2(2) dim1(3)]);
S_A = kron(kr(C1,B2),eye(dim2(1)));
S_B = J_12'*kron(kr(C1,A2),eye(dim2(2)));
S_C = J_13'*kron(kr(B2,A2),eye(dim1(3)));
ind = [1:dim2(1):dim2(1)*F dim2(1)*F+1:dim2(2):(dim2(1)+dim2(2))*F];
M1 = eye((dim2(1)+dim2(2)+dim1(3))*F); M1(ind,:) = []; %For omega
Hmat = [S_A S_B S_C]*M1'; Hmat = sparse(Hmat);


for q=1:2:36
    
    phi = gauss_kernel(q);
    H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);
        
    for d=1:6
        
        clear FIM1 FIM2 FIM CCRB
    
        S = eye(dim2(1)); S = S(1:d:end,:);
        P1 = S*H; P1 = P1(1:dim1(1),:); P2 = P1;

        A1 = P1*A2; B1 = P2*B2; C2 = Pm*C1;
        X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2}); %Groundtruth tensors
        
        for s=1:length(sigma_n2)
            [~,FIM2(:,:,s)] = crb_uncoupled(A2,B2,C2,sigma_n2(s));
        end
        for s=1:length(sigma_n1)
            [~, FIM1(:,:,s)] = crb_uncoupled2(A1,B1,C1,sigma_n1(s));
            FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2); %For CCRB
        end
        
        H3 = [-kron(eye(F),Pm) zeros(dim2(3)*F,(dim1(1)+dim1(2))*F) eye(dim2(3)*F) zeros(dim2(3)*F,(dim2(1)+dim2(2)-2)*F); ...
        zeros((dim1(1)+dim1(2))*F,dim1(3)*F) eye((dim1(1)+dim1(2))*F) zeros((dim1(1)+dim1(2))*F,dim2(3)*F) -blkdiag(kron(eye(F),P1(:,2:end)),kron(eye(F),P2(:,2:end)))];
        U = null(H3);
       
        
        for s=1:length(sigma_n1)
            
            [q d SNR1(s)]
            CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
            
            CCRB4(:,:,s) = [CCRB((sum(dim1)-2)*F+1:end,(sum(dim1)-2)*F+1:end,s) CCRB((sum(dim1)-2)*F+1:end,1:(sum(dim1)-2)*F,s);
                    CCRB(1:(sum(dim1)-2)*F,(sum(dim1)-2)*F+1:end,s) CCRB(1:(sum(dim1)-2)*F,1:(sum(dim1)-2)*F,s)];
            temp1(:,:,s) = CCRB4(dim2(3)*F+1:(sum(dim2)+dim1(3)-2)*F,dim2(3)*F+1:(sum(dim2)+dim1(3)-2)*F,s);
            CCRB_rec(q,d,s) = sum(diag(sparse(Hmat)*temp1(:,:,s)*sparse(Hmat')));

        end
    end
end

%% Process

SNR = 5:5:60; SNR1 = SNR;

load('allbounds_rank3.mat')
CCRB_rec = CCRB_rec(1:2:end,:,:);
for s=1:size(CCRB_rec,3)
    minval = min(CCRB_rec(:,:,s),[],'all');
    [row_q1(s),col_d1(s)] = find(CCRB_rec(:,:,s) == minval);
    bound1(s) = CCRB_rec(row_q1(s),col_d1(s),s);
end

load('allbounds_rank16.mat')
CCRB_rec = CCRB_rec(1:2:end,:,:);
for s=1:size(CCRB_rec,3)
    minval = min(CCRB_rec(:,:,s),[],'all');
    [row_q2(s),col_d2(s)] = find(CCRB_rec(:,:,s) == minval);
    bound2(s) = CCRB_rec(row_q2(s),col_d2(s),s);
end

figure(1)
subplot(1,2,1)
plot(SNR1,col_d1,'k+'); hold on; plot(SNR1,col_d2,'ro');
title('Optimal $d$','Interpreter','latex')
ylabel('Downsampling ratio $d$','Interpreter','latex')
xlabel('$\mathbf{SNR}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
xlim([5 60])
ylim([0 7])
subplot(1,2,2)
plot(SNR1,row_q1,'k+'); hold on; plot(SNR1,row_q2,'ro');
title('Optimal $q$','Interpreter','latex')
ylabel('Filter size $q$','Interpreter','latex')
xlabel('$\mathbf{SNR}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
xlim([5 60])
ylim([0 20])
saveas(gcf,'figures/fig_optimal_qd.fig')

figure(2)
subplot(1,2,1)
semilogy(SNR,bound1,'k--','LineWidth',1)
title('$N=3$','Interpreter','latex')
ylabel('$\mathbf{CCRB(\psi)}$','Interpreter','latex')
xlabel('$\mathbf{SNR}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
xlim([5 60])
subplot(1,2,2)
semilogy(SNR,bound2,'k--','LineWidth',1)
title('$N=16$','Interpreter','latex')
ylabel('$\mathbf{CCRB(\psi)}$','Interpreter','latex')
xlabel('$\mathbf{SNR}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
xlim([5 60])
saveas(gcf,'figures/ccrb_optimal_qd.fig')
