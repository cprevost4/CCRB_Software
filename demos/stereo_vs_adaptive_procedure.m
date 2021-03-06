clear all
close all
clc
rng(107)

prompt = "Select the rank (the preprint uses ranks 10,12,14)";
num = input(prompt);
eval(sprintf("F=%d",num));
fprintf("Please be patient, this code is quite long due to the dimensions.")


%% Groundtruth data

dim1 = [4,4,20]; dim2 = [16,16,10]; %Dimensions oF CP model

SNR1 = 5:5:60; %Noise on first tensor
for s=1:length(SNR1)
    sigma_n1(s) = 10^(-SNR1(s)/10);
end
SNR2 = 40; sigma_n2 = 10^(-SNR2/10); %Noise on second tensor

q = 3; phi = gauss_kernel(q); 
H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);
d = 4; S = eye(d*dim1(1)); S = S(1:d:end,:);
P1 = S*H; P2 = P1;
Pm = eye(dim1(3));
Pm(1:end-1,2:end) = Pm(1:end-1,2:end) + eye(dim1(3)-1); Pm = Pm(1:2:end,:);
Pm = Pm/2;

A2 = randn(dim2(1),F); B2 = randn(dim2(2),F); C1 = randn(dim1(3),F);
A2(1,:) = 1; B2(1,:) = 1;
alpha = P1(1,:)*A2; beta = P2(1,:)*B2;
A1 = P1*A2*diag(1./alpha);
B1 = P2*B2*diag(1./beta);
C2 = Pm*C1*diag(1./(alpha.*beta));

X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2}); %Groundtruth tensors

%% CRB

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
    FIM(:,:,s) = blkdiag(FIM1(:,:,s),FIM2); %For CCRB
end

%% CCRB 

[~, H3] = transformation_jac(A2,B2,C1,P1,P2,Pm,alpha,beta,dim1,dim2);
U = null(H3);

for s=1:length(sigma_n1)
    CCRB(:,:,s) = U*inv(U'*FIM(:,:,s)*U)*U';
    CCRB2_C1(s) = sum(diag(CCRB(1:dim1(3)*F,1:dim1(3)*F,s)));
    CCRB2_AB1(s) = sum(diag(CCRB(dim1(3)*F+1:(sum(dim1)-2)*F,dim1(3)*F+1:(sum(dim1)-2)*F,s)));
    CCRB2_C2(s) = sum(diag(CCRB((sum(dim1)-2)*F+1:(dim2(3)+sum(dim1)-2)*F,(sum(dim1)-2)*F+1:(dim2(3)+sum(dim1)-2)*F,s)));
    CCRB2_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*F+1:end,(dim2(3)+sum(dim1)-2)*F+1:end,s))); 
    CCRB1_psi1(s) = CCRB2_AB2(s)+CCRB2_C1(s);
    CCRB1_psi2(s) = CCRB2_AB1(s)+CCRB2_C2(s);
end

%% PRE-ALLOCATION FOR SPEED

Nreal = 1; Ninit = 50; Niter_c = 5000; Niter = 1;

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

%% STEREO

P = Nreal*Ninit*length(sigma_n1); p=0; %Progress indicator 

%c = logspace(0,log10(sigma_n2),5);

for s=1:length(sigma_n1)
    
    c = logspace(0,log10(sigma_n1(s)),5);
    
    %lambda = 1; %c = sigma_n1(s);
    %lambda = (sigma_n2/sigma_n1(s))^2; %true
    %lambda = (sigma_n1(s)/sigma_n2)^2; % inverse
     
    for n=1:Nreal
        
        %Noisy tensors
        Y1 = X1+sigma_n1(s)*randn(dim1);
        Y2 = X2+sigma_n2*randn(dim2);
        %Initial objective values
        obj = 10^50; obj1 = 10^50; obj2 = 10^50; obj3 = 10^50;
        
        for i=1:Ninit
            
            %Initialization
            U2_0 = cpd(Y2,F);
            A2_0 = U2_0{1}; B2_0 = U2_0{2}; C2_0 = U2_0{3};
            A2_0(1,:) = 1; B2_0(1,:) = 1;
            temp = kr(P1*A2_0,P2*B2_0);
            C1_0=(temp\tens2mat(Y1,[],3))';
% U1_0 = cpd(Y1,F); U2_0 = cpd(Y2,F);
%             A1_0 = U1_0{1}; B1_0 = U1_0{2}; C1_0 = U1_0{3};
%             A2_0 = U2_0{1}; B2_0 = U2_0{2}; C2_0 = U2_0{3};
%       
%             C1_0(1,:) = 1; 
%             A1_0(1,:) = 1; B1_0(1,:) = 1;
%             C2_0(1,:) = 1; 
%             A2_0(1,:) = 1; B2_0(1,:) = 1;
            
%             A2_0 = randn(size(A2)); B2_0 = randn(size(B2)); C2_0 = randn(size(C2));
%             C1_0 = randn(size(C1));
%             A2_0(1,:) = 1; B2_0(1,:) = 1;
            
            %for l=1:length(c)
                lambda = (sigma_n2/sigma_n1(s))^2; %true
                %lambda(s,l) = (sigma_n2/c(l))^2;
                % ALS algorithm, coupled
                [A2_n,B2_n,C1_n,cost] = STEREO( Y1,Y2,P1,P2,Pm,Niter_c,lambda,A2_0,B2_0,C1_0,C2_0, [sigma_n1(s) sigma_n2]);
                %A2_0 = A2_n; B2_0 = B2_n; C1_0 = C1_n; C2_0 = Pm*C1_n;
           % end
            
            %Pick best init
            if cost(end)<obj
                obj = cost(end);
                C1_c = C1_n; A2_c = A2_n; B2_c = B2_n;
            end
            
            p=p+1; clc
            fprintf('Progress %g %%',(p/P)*100)
        end
           
        %Permute result of STEREO
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
        
        %Squared errors
        se_C1_c(:,n,s) = (C1(:)-C1_c(:)).^2;
        err = [A2(2:end,:) - A2_c(2:end,:) B2(2:end,:) - B2_c(2:end,:)];
        se_AB2_c(:,n,s) = (err(:)).^2; 
        se_C2_c(:,n,s) = (C2(:)-C2_c(:)).^2;
        err = [A1(2:end,:) - A1_c(2:end,:) B1(2:end,:) - B1_c(2:end,:)];
        se_AB1_c(:,n,s) = (err(:)).^2;  
    end  
    
    %MSE
    
    mse_C1_c(s) = sum(mean(se_C1_c(:,:,s),2));
    mse_AB2_c(s) = sum(mean(se_AB2_c(:,:,s),2));
    mse_C2_c(s) = sum(mean(se_C2_c(:,:,s),2));
    mse_AB1_c(s) = sum(mean(se_AB1_c(:,:,s),2));

    mse_psi1_c(s) = mse_AB2_c(s)+mse_C1_c(s);
    mse_psi2_c(s) = mse_AB1_c(s)+mse_C2_c(s);


end

%% Adaptive procedure

P = Nreal*Ninit*length(sigma_n1); p=0; %Progress indicator 

%c = logspace(0,log10(sigma_n2),5);

for s=1:length(sigma_n1)
    
    c = logspace(0,log10(sigma_n1(s)),5);
    
    %lambda = 1; %c = sigma_n1(s);
    %lambda = (sigma_n2/sigma_n1(s))^2; %true
    %lambda = (sigma_n1(s)/sigma_n2)^2; % inverse
     
    for n=1:Nreal
        
        %Noisy tensors
        Y1 = X1+sigma_n1(s)*randn(dim1);
        Y2 = X2+sigma_n2*randn(dim2);
        %Initial objective values
        obj = 10^50; obj1 = 10^50; obj2 = 10^50; obj3 = 10^50;
        
        for i=1:Ninit
            
            %Initialization
            U2_0 = cpd(Y2,F);
            A2_0 = U2_0{1}; B2_0 = U2_0{2}; C2_0 = U2_0{3};
            A2_0(1,:) = 1; B2_0(1,:) = 1;
            temp = kr(P1*A2_0,P2*B2_0);
            C1_0=(temp\tens2mat(Y1,[],3))';
% U1_0 = cpd(Y1,F); U2_0 = cpd(Y2,F);
%             A1_0 = U1_0{1}; B1_0 = U1_0{2}; C1_0 = U1_0{3};
%             A2_0 = U2_0{1}; B2_0 = U2_0{2}; C2_0 = U2_0{3};
%       
%             C1_0(1,:) = 1; 
%             A1_0(1,:) = 1; B1_0(1,:) = 1;
%             C2_0(1,:) = 1; 
%             A2_0(1,:) = 1; B2_0(1,:) = 1;
            
%             A2_0 = randn(size(A2)); B2_0 = randn(size(B2)); C2_0 = randn(size(C2));
%             C1_0 = randn(size(C1));
%             A2_0(1,:) = 1; B2_0(1,:) = 1;
            
            for l=1:length(c)
                %lambda = (sigma_n2/sigma_n1(s))^2; %true
                lambda(s,l) = (sigma_n2/c(l))^2;
                % ALS algorithm, coupled
                [A2_n,B2_n,C1_n,cost] = STEREO( Y1,Y2,P1,P2,Pm,Niter_c/5,lambda(s,l),A2_0,B2_0,C1_0,C2_0, [sigma_n1(s) sigma_n2]);
                A2_0 = A2_n; B2_0 = B2_n; C1_0 = C1_n; C2_0 = Pm*C1_n;
            end
            
            %Pick best init
            if cost(end)<obj
                obj = cost(end);
                C1_c = C1_n; A2_c = A2_n; B2_c = B2_n;
            end
            
            p=p+1; clc
            fprintf('Progress %g %%',(p/P)*100)
        end
           
        %Permute result of STEREO
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
        
        %Squared errors
        se_C1_c(:,n,s) = (C1(:)-C1_c(:)).^2;
        err = [A2(2:end,:) - A2_c(2:end,:) B2(2:end,:) - B2_c(2:end,:)];
        se_AB2_c(:,n,s) = (err(:)).^2; 
        se_C2_c(:,n,s) = (C2(:)-C2_c(:)).^2;
        err = [A1(2:end,:) - A1_c(2:end,:) B1(2:end,:) - B1_c(2:end,:)];
        se_AB1_c(:,n,s) = (err(:)).^2;  
    end  
    
    %MSE
    
    mse_C1_c(s) = sum(mean(se_C1_c(:,:,s),2));
    mse_AB2_c(s) = sum(mean(se_AB2_c(:,:,s),2));
    mse_C2_c(s) = sum(mean(se_C2_c(:,:,s),2));
    mse_AB1_c(s) = sum(mean(se_AB1_c(:,:,s),2));

    mse_psi1_c2(s) = mse_AB2_c(s)+mse_C1_c(s);
    mse_psi2_c2(s) = mse_AB1_c(s)+mse_C2_c(s);


end


%% FIGURES

figure(1)
semilogy(SNR1,CCRB1_psi1+CCRB1_psi2,'k--','LineWidth',1); hold on
semilogy(SNR1,mse_psi1_c+mse_psi2_c,'+','MarkerSize',10); hold on
semilogy(SNR1,mse_psi1_c2+mse_psi2_c2,'rd'); hold on
legend('CCRB','$\lambda = \frac{\sigma_1^2}{\sigma_2^2}$','Adaptive $\lambda$','Interpreter','latex')
title(sprintf('Performance for $N=%d$',F),'Interpreter','latex')
xlabel('$\mathbf{SNR}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',14)
xlim([5 60])
saveas(gcf,'figures/fig_adaptive_procedure.fig')

