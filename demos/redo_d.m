clear all
close all
clc
rng(107)

prompt = "Select the rank (between 3 and 16)";
num = input(prompt);
eval(sprintf("F=%d",num));
fprintf("Please be patient, this code is quite long due to the dimensions.")

%% Run simulations

dim1 = [6,6,16]; dim2 = [36,36,8];

SNR1 = [5 15 25 30 45 60];
for s=1:length(SNR1)
    sigma_n1(s) = F*10^(-SNR1(s)/10);
end
SNR2 = 20; sigma_n2 = F*10^(-SNR2/10); 

Pm = eye(dim1(3));
Pm(1:end-1,2:end) = Pm(1:end-1,2:end) + eye(dim1(3)-1); Pm = Pm(1:2:end,:);
Pm = Pm/2;
q=3; phi = gauss_kernel(q);
H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);

A2 = randn(dim2(1),F); B2 = randn(dim2(2),F); C1 = randn(dim1(3),F);
A2(1,:) = 1; B2(1,:) = 1;

[J_12,J_13,~] = vec_unfold([dim2(1)-1 dim2(2)-1 dim1(3)]);
S_A = kron(kr(C1,B2(2:end,:)),eye(dim2(1)-1));
S_B = J_12'*kron(kr(C1,A2(2:end,:)),eye(dim2(2)-1));
S_C = J_13'*kron(kr(B2(2:end,:),A2(2:end,:)),eye(dim1(3)));
Hmat = [S_A S_B S_C]; Hmat = sparse(Hmat);

for d=1:6
    S = eye(dim2(1)); S = S(1:d:end,:);
    P1 = S*H; P1 = P1(1:dim1(1),:); P2 = P1;
    cond_P1(d) = cond(P1);
     
    A1 = P1*A2; B1 = P2*B2; C2 = Pm*C1;
    X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2}); 

    for s=1:length(sigma_n2)
        [~,FIM2(:,:,s)] = crb_uncoupled(A2,B2,C2,sigma_n2(s));
    end
    for s=1:length(sigma_n1)
        [~, FIM1(:,:,s)] = crb_uncoupled2(A1,B1,C1,sigma_n1(s));
        FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2); 
    end

    H3 = [-kron(eye(F),Pm) zeros(dim2(3)*F,(dim1(1)+dim1(2))*F) eye(dim2(3)*F) zeros(dim2(3)*F,(dim2(1)+dim2(2)-2)*F); ...
    zeros((dim1(1)+dim1(2))*F,dim1(3)*F) eye((dim1(1)+dim1(2))*F) zeros((dim1(1)+dim1(2))*F,dim2(3)*F) -blkdiag(kron(eye(F),P1(:,2:end)),kron(eye(F),P2(:,2:end)))];
    U = null(H3);
    for s=1:length(sigma_n1)
        CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
        CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F,1:dim1(3)*F,s)));
        CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1))*F+1:end,(dim2(3)+sum(dim1))*F+1:end,s)));
        CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);

        CCRB4(:,:,s) = [CCRB((sum(dim1)-2)*F+1:end,(sum(dim1)-2)*F+1:end,s) CCRB((sum(dim1)-2)*F+1:end,1:(sum(dim1)-2)*F,s);
                CCRB(1:(sum(dim1)-2)*F,(sum(dim1)-2)*F+1:end,s) CCRB(1:(sum(dim1)-2)*F,1:(sum(dim1)-2)*F,s)];
        temp1(:,:,s) = CCRB4(dim2(3)*F+1:(sum(dim2)+dim1(3)-2)*F,dim2(3)*F+1:(sum(dim2)+dim1(3)-2)*F,s);
        CCRB_rec(d,s) = sum(diag(sparse(Hmat)*temp1(:,:,s)*sparse(Hmat'))); 
    end

    bound = CCRB_rec;
        
end
 

figure(1)
sgtitle('Performance bounds for $N=16$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
for s=1:3
    subplot(1,2,1)
    semilogy(1:6,bound(:,s),'--','LineWidth',1); hold on
    xlim([1 6])
    ylabel('$\mathbf{CCRB(\psi)}$','Interpreter','latex')
    xlabel('Downsampling ratio $d$','Interpreter','latex')
    set(gca,'FontName','Times','FontSize',16)
end
for s=4:6
    subplot(1,2,2)
    semilogy(1:6,bound(:,s),'--','LineWidth',1); hold on
    xlim([1 6])
    ylabel('$\mathbf{CCRB(\psi)}$','Interpreter','latex')
    xlabel('Downsampling ratio $d$','Interpreter','latex')
    set(gca,'FontName','Times','FontSize',16)
end
saveas(gcf,'figures/fig_impact_downsamplingration.fig')


