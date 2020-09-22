clear all
close all
clc
rng(107) %For reproducibility

%% Groundtruth data

%Dimensions of CP model
dim1 = [6,6,16]; dim2 = [18,18,8]; F = 3; 

%Noise level
SNR1 = 5:5:60; 
for s=1:length(SNR1)
    sigma_n1(s) = F*10^(-SNR1(s)/10);
end
SNR2 = 20; sigma_n2 = F*10^(-SNR2/10);

%Generate degradation matrices
q = 3; phi = gauss_kernel(q);
H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);
d = 3; S = eye(d*dim1(1)); S = S(1:d:end,:);
P1 = S*H; P2 = P1;
Pm = eye(dim1(3));
Pm(1:end-1,2:end) = Pm(1:end-1,2:end) + eye(dim1(3)-1); Pm = Pm(1:2:end,:);
Pm = Pm/2;

%Generate low-rank factors
A2 = randn(dim2(1),F); B2 = randn(dim2(2),F); C1 = randn(dim1(3),F);
A2(1,:) = 1; B2(1,:) = 1;
alpha = P1(1,:)*A2; beta = P2(1,:)*B2;
A1 = P1*A2*diag(1./alpha);
B1 = P2*B2*diag(1./beta);
C2 = Pm*C1*diag(1./(alpha.*beta));

%Groundtruth tensors
X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2});

%% Uncoupled CRB

%CRB for Y2
for s=1:length(sigma_n2)
    [CRB2(:,:,s),FIM2(:,:,s)] = crb_uncoupled(A2,B2,C2,sigma_n2(s));
    CRB_C2(s) = sum(diag(CRB2(1:dim2(3)*F,1:dim2(3)*F,s)));
    CRB_AB2(s) = sum(diag(CRB2(dim2(3)*F+1:end,dim2(3)*F+1:end,s)));
end

%CRB for Y1
for s=1:length(sigma_n1)
    [CRB1(:,:,s), FIM1(:,:,s)] = crb_uncoupled(A1,B1,C1,sigma_n1(s));
    CRB_C1(s) = sum(diag(CRB1(1:dim1(3)*F,1:dim1(3)*F,s)));
    CRB_AB1(s) = sum(diag(CRB1(dim1(3)*F+1:end,dim1(3)*F+1:end,s)));

    CRB_psi1(s) = CRB_AB2+CRB_C1(s);
    CRB_psi2(s) = CRB_AB1(s)+CRB_C2;

    FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2);
end

%% Fully coupled CCRB

%Compute jacobian of the (non-linear) constraints on the parameters
[~, H3] = transformation_jac(A2,B2,C1,P1,P2,Pm,alpha,beta,dim1,dim2);

%Basis for null-space of H3
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

%% Blind-CCRB

%Jacobian of the constraints
dC2_dC1 = kron(diag(1./(alpha.*beta)),Pm);
dC2_dA2 = -kron(diag(1./(beta.*(alpha.^2))),Pm)*kr(eye(F),C1)*kron(eye(F),P1(1,2:end));
dC2_dB2 = -kron(diag(1./(alpha.*(beta.^2))),Pm)*kr(eye(F),C1)*kron(eye(F),P2(1,2:end));
H3 = [-dC2_dC1 zeros(dim2(3)*F,(dim1(1)+dim1(2)-2)*F) eye(dim2(3)*F) -dC2_dA2 -dC2_dB2];
%Basis for null-space of H3
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

%% Pre-allocation (for speed purposes)

%Specify the values for simulations
Nreal = 200; Ninit = 10; Niter_c = 5000; Niter = 1000;

%Pre-allocation for MSE
se_C1 = zeros(dim1(3)*F,Nreal,length(SNR1));
se_AB2 = zeros((dim2(1)+dim2(2)-2)*F,Nreal,length(SNR1));
se_C2 = zeros(dim2(3)*F,Nreal,length(SNR1));
se_AB1 = zeros((dim1(1)+dim1(2)-2)*F,Nreal,length(SNR1));
se_C1_c = zeros(dim1(3)*F,Nreal,length(SNR1));
se_AB2_c = zeros((dim2(1)+dim2(2)-2)*F,Nreal,length(SNR1));
se_C2_c = zeros(dim2(3)*F,Nreal,length(SNR1));
se_AB1_c = zeros((dim1(1)+dim1(2)-2)*F,Nreal,length(SNR1));
mse_C1 = zeros(1,length(SNR1));
mse_AB2 = zeros(1,length(SNR1));
mse_C2 = zeros(1,length(SNR1));
mse_AB1 = zeros(1,length(SNR1));
mse_C1_c = zeros(1,length(SNR1));
mse_AB2_c = zeros(1,length(SNR1));
mse_C2_c = zeros(1,length(SNR1));
mse_AB1_c = zeros(1,length(SNR1));
mse_psi1 = zeros(1,length(SNR1));
mse_psi2 = zeros(1,length(SNR1));
mse_psi1_c = zeros(1,length(SNR1));
mse_psi2_c = zeros(1,length(SNR1));

%% Simulations

%Progress indicator 
P = Nreal*Ninit*length(sigma_n1); p=0; 

for s=1:length(sigma_n1)
    for n=1:Nreal
        
        %Generate noisy tensors
        Y1 = X1+sigma_n1(s)*randn(dim1);
        Y2 = X2+sigma_n2*randn(dim2);
        %Initial objective values
        obj = 10^50; obj1 = 10^50; obj2 = 10^50; obj3 = 10^50;
        
        for i=1:Ninit
            
            %Initialization
            U1_0 = cpd(Y1,F); U2_0 = cpd(Y2,F);
            A1_0 = U1_0{1}; B1_0 = U1_0{2}; C1_0 = U1_0{3};
            A2_0 = U2_0{1}; B2_0 = U2_0{2}; C2_0 = U2_0{3};
      
            C1_0(1,:) = 1; 
            A1_0(1,:) = 1; B1_0(1,:) = 1;
            C2_0(1,:) = 1; 
            A2_0(1,:) = 1; B2_0(1,:) = 1;
            
            %ALS algorithm, uncoupled
            [A1_nu,B1_nu,C1_nu,cost1] = als(Y1,Niter,A1_0,B1_0,C1_0);
            [A2_nu,B2_nu,C2_nu,cost2] = als(Y2,Niter,A2_0,B2_0,C2_0);
            
            %Pick best init
            if cost1(end)<obj1
                obj1 = cost1(end);
                A1_u = A1_nu; B1_u = B1_nu; C1_u = C1_nu;
            end
            if cost2(end)<obj2
                obj2 = cost2(end);
                A2_u = A2_nu; B2_u = B2_nu; C2_u = C2_nu;
             end
            
            %-----------------------
            
            % ALS algorithm, coupled
            [A2_n,B2_n,C1_n,cost] = STEREO( Y1,Y2,P1,P2,Pm,Niter_c,(sigma_n2/sigma_n1(s))^2,A2_u,B2_u,C1_u,C2_u, [sigma_n1(s) sigma_n2]);
            
            %Pick best init
            if cost(end)<obj
                obj = cost(end);
                C1_c = C1_n; A2_c = A2_n; B2_c = B2_n;
            end

            % Blind coupled ALS algo
            [A2_n,B2_n,C1_n,A1_n,B1_n,C2_n,cost] = Blind_STEREO(Y1,Y2,Pm,Niter_c,A2_u,B2_u,A1_u,B1_u,Pm*C1_u,(sigma_n2/sigma_n1(s))^2,[sigma_n1(s) sigma_n2]);
            
            %Pick best init
            if cost(end)<obj3
                obj = cost(end);
                C1_b = C1_n; A2_b = A2_n; B2_b = B2_n; A1_b = A1_n; B1_b = B1_n;
                C2_b = C2_n;
            end
            
            %-----------------------

            p=p+1; clc
            fprintf('Progress %g %%',(p/P)*100)
        end

        
        
        %Correct ambiguity - uncoupled ALS
        mat = C1_u'*C1; ind = [];
        tmp = zeros(size(A1)); tmp2 = zeros(size(B1)); tmp3 = zeros(size(C1));
        for i=1:F
              ind(i) = find(mat(:,i)==max(mat(:,i)));
              tmp(:,i) = A1_u(:,ind(i)); tmp2(:,i) = B1_u(:,ind(i));
              tmp3(:,i) = C1_u(:,ind(i));
        end
        A1_u = tmp; B1_u = tmp2; C1_u = tmp3;
        mat = A2_u'*A2; ind = [];
        tmp = zeros(size(A2)); tmp2 = zeros(size(B2)); tmp3 = zeros(size(C2));
        for i=1:F
              ind(i) = find(mat(:,i)==max(mat(:,i)));
              tmp(:,i) = A2_u(:,ind(i)); tmp2(:,i) = B2_u(:,ind(i));
              tmp3(:,i) = C2_u(:,ind(i));
        end
        A2_u = tmp; B2_u = tmp2; C2_u = tmp3;
        Y_hat_u = cpdgen({A2_u(2:end,:),B2_u(2:end,:),C1_u});
        
        %-----------------------
           
        %Correct ambiguity - STEREO
        mat = A2_c'*A2; ind = [];
        tmp = zeros(size(A2)); tmp2 = zeros(size(B2)); tmp3 = zeros(size(C1));
        for i=1:F
              ind(i) = find(mat(:,i)==max(mat(:,i)));
              tmp(:,i) = A2_c(:,ind(i)); tmp2(:,i) = B2_c(:,ind(i));
              tmp3(:,i) = C1_c(:,ind(i));
        end
        A2_c = tmp; B2_c = tmp2; C1_c = tmp3;
        
        %Generate degraded matrices - STEREO
        A1_c =     (P1*A2_c).*repmat(1./alpha,dim1(1),1);
        B1_c =     (P2*B2_c).*repmat(1./beta,dim1(2),1);
        C2_c =     (Pm*C1_c).*repmat(1./(alpha.*beta),dim2(3),1);
        
        C2_c = C2_c.*repmat(1./(A1_c(1,:).*B1_c(1,:)),dim2(3),1);
        A1_c = A1_c.*repmat(1./A1_c(1,:),dim1(1),1);
        B1_c = B1_c.*repmat(1./B1_c(1,:),dim1(2),1);
        
        %-----------------------

        %Correct ambiguity - Blind-STEREO
        mat = A1_b'*A1; ind = [];
        tmp = zeros(size(A1)); tmp2 = zeros(size(B1));
        for i=1:F
              ind(i) = find(mat(:,i)==max(mat(:,i)));
              tmp(:,i) = A1_b(:,ind(i)); tmp2(:,i) = B1_b(:,ind(i));
              tmp3(:,i) = C1_b(:,ind(i));
        end
        A1_b = tmp; B1_b = tmp2; C1_b = tmp3; 
        mat = A2_b'*A2; ind = [];
        tmp = zeros(size(A2)); tmp2 = zeros(size(B2)); tmp3 = zeros(size(C1)); tmp4 = zeros(size(C2));
        for i=1:F
              ind(i) = find(mat(:,i)==max(mat(:,i)));
              tmp(:,i) = A2_b(:,ind(i)); tmp2(:,i) = B2_b(:,ind(i));
              tmp4(:,i) = C2_b(:,ind(i)); 
        end
        A2_b = tmp; B2_b = tmp2;
        C2_b = tmp4;
        
        %-----------------------
        
        %Squared errors
        se_C1(:,n,s) = (C1(:)-C1_u(:)).^2;
        err = [A2(2:end,:) - A2_u(2:end,:) B2(2:end,:) - B2_u(2:end,:)];
        se_AB2(:,n,s) = (err(:)).^2; 
        se_C2(:,n,s) = (C2(:)-C2_u(:)).^2;
        err = [A1(2:end,:) - A1_u(2:end,:) B1(2:end,:) - B1_u(2:end,:)];
        se_AB1(:,n,s) = (err(:)).^2;    
        
        se_C1_c(:,n,s) = (C1(:)-C1_c(:)).^2;
        err = [A2(2:end,:) - A2_c(2:end,:) B2(2:end,:) - B2_c(2:end,:)];
        se_AB2_c(:,n,s) = (err(:)).^2; 
        se_C2_c(:,n,s) = (C2(:)-C2_c(:)).^2;
        err = [A1(2:end,:) - A1_c(2:end,:) B1(2:end,:) - B1_c(2:end,:)];
        se_AB1_c(:,n,s) = (err(:)).^2;  

        se_C1_b(:,n,s) = (C1(:)-C1_b(:)).^2;
        err = [A2(2:end,:) - A2_b(2:end,:) B2(2:end,:) - B2_b(2:end,:)];
        se_AB2_b(:,n,s) = (err(:)).^2; 
        se_C2_b(:,n,s) = (C2(:)-C2_b(:)).^2;
        err = [A1(2:end,:) - A1_b(2:end,:) B1(2:end,:) - B1_b(2:end,:)];
        se_AB1_b(:,n,s) = (err(:)).^2; 

    end  
    
    %MSE
    mse_C1(s) = sum(mean(se_C1(:,:,s),2));
    mse_AB2(s) = sum(mean(se_AB2(:,:,s),2));
    mse_C2(s) = sum(mean(se_C2(:,:,s),2));
    mse_AB1(s) = sum(mean(se_AB1(:,:,s),2));
    
    mse_C1_c(s) = sum(mean(se_C1_c(:,:,s),2));
    mse_AB2_c(s) = sum(mean(se_AB2_c(:,:,s),2));
    mse_C2_c(s) = sum(mean(se_C2_c(:,:,s),2));
    mse_AB1_c(s) = sum(mean(se_AB1_c(:,:,s),2));

    mse_C1_b(s) = sum(mean(se_C1_b(:,:,s),2));
    mse_AB2_b(s) = sum(mean(se_AB2_b(:,:,s),2));
    mse_C2_b(s) = sum(mean(se_C2_b(:,:,s),2));
    mse_AB1_b(s) = sum(mean(se_AB1_b(:,:,s),2));
    
    mse_psi1(s) = mse_AB2(s)+mse_C1(s);
    mse_psi2(s) = mse_AB1(s)+mse_C2(s);
    mse_psi1_c(s) = mse_AB2_c(s)+mse_C1_c(s);
    mse_psi2_c(s) = mse_AB1_c(s)+mse_C2_c(s);
    mse_psi1_b(s) = mse_AB2_b(s)+mse_C1_b(s);
    mse_psi2_b(s) = mse_AB1_b(s)+mse_C2_b(s);

end

%% Figure

figure(1)
subplot(1,2,1)
semilogy(SNR1,CRB_psi1,'k-','LineWidth',1); hold on
semilogy(SNR1,mse_psi1,'k*','MarkerSize',9); hold on
semilogy(SNR1,CCRB1_psi1,'b--','LineWidth',1); hold on;
semilogy(SNR1,mse_psi1_c,'bo','MarkerSize',9); hold on
semilogy(SNR1,CCRB2_psi1,'r:','LineWidth',2); hold on;
semilogy(SNR1,mse_psi1_b,'rd','MarkerSize',9); hold on
title('Performance bounds for $\widetilde{\mathbf{\theta}}$','Interpreter','latex')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
legend('Uncoupled CRB','MSE - Uncoupled ALS', 'Fully-coupled CCRB', 'MSE - STEREO', 'Blind-CCRB', 'MSE - Blind-STEREO')
set(gca,'FontName','Times','FontSize',16)
xlim([5 60])
    
subplot(1,2,2)
semilogy(SNR1,CRB_psi2,'k-','LineWidth',1); hold on
semilogy(SNR1,mse_psi2,'k*','MarkerSize',9); hold on
semilogy(SNR1,CCRB1_psi2,'b--','LineWidth',1); hold on
semilogy(SNR1,mse_psi2_c,'bo','MarkerSize',9); hold on
semilogy(SNR1,CCRB2_psi2,'r:','LineWidth',2); hold on;
semilogy(SNR1,mse_psi2_b,'rd','MarkerSize',9); hold on
title('Performance bounds for $\widetilde{\mathbf{\phi}}$','Interpreter','latex')
xlabel('SNR on $\mathcal{Y}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
legend('Uncoupled CRB','MSE - Uncoupled ALS', 'Fully-coupled CCRB', 'MSE - STEREO', 'Blind-CCRB', 'MSE - Blind-STEREO')
xlim([5 60])
