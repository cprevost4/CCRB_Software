clear all
close all
clc
rng(107) %For reproducibility

%% Groundtruth data

%Dimensions oF CP model
dim1 = [10,10,15]; dim2 = [15,15,9]; F = 3;

%Noise on tensors
SNR1 = 5:5:40; 
for s=1:length(SNR1)
    sigma_n1(s) = F*10^(-SNR1(s)/10);
end
SNR2 = 20; sigma_n2 = F*10^(-SNR2/10)*ones(1,length(SNR1)); 

%Degradation matrices
P1 = eye(dim2(1)); P1 = P1(1:dim1(1),:); 
P2 = eye(dim2(2)); P2 = P2(1:dim1(2),:); 
Pm = eye(dim1(3)); Pm = Pm(1:dim2(3),:);

%CP factors
A2 = randn(dim2(1),F); B2 = randn(dim2(2),F); C1 = randn(dim1(3),F);
A1 = P1*A2; B1 = P2*B2; C2 = Pm*C1;
A1(1,:) = ones(1,F); B1(1,:) = ones(1,F);
A2(1,:) = ones(1,F); B2(1,:) = ones(1,F);

%Groundtruth tensors
X1 = cpdgen({A1,B1,C1}); X2 = cpdgen({A2,B2,C2});

%% Uncoupled CRB

%CRB for X1
for s=1:length(sigma_n1)
    CRB1(:,:,s) = crb_uncoupled(A1,B1,C1,sigma_n1(s));
    CRB_C1(s) = sum(diag(CRB1(1:dim1(3)*F,1:dim1(3)*F,s)));
    CRB_AB1(s) = sum(diag(CRB1(dim1(3)*F+1:end,dim1(3)*F+1:end,s)));
end
%CRB for X2
for s=1:length(sigma_n2)
    CRB2(:,:,s) = crb_uncoupled(A2,B2,C2,sigma_n2(s));
    CRB_C2(s) = sum(diag(CRB2(1:dim2(3)*F,1:dim2(3)*F,s)));
    CRB_AB2(s) = sum(diag(CRB2(dim2(3)*F+1:end,dim2(3)*F+1:end,s)));
end


%% CCRB

%Jacobian of the constraints
P = [kron(eye(F),P1) zeros(dim1(1)*F,dim2(2)*F);
    zeros(dim1(2)*F,dim2(1)*F) kron(eye(F),P2)];
ind = [1:dim1(1):dim1(1)*F dim1(1)*F+1:dim1(2):(dim1(1)+dim1(2))*F];
M1 = eye((dim1(1)+dim1(2))*F); M1(ind,:) = []; 
ind = [1:dim2(1):dim2(1)*F dim2(1)*F+1:dim2(2):(dim2(1)+dim2(2))*F];
M2 = eye((dim2(1)+dim2(2))*F); M2(ind,:) = []; 
P_tilde = M1*P*M2';
G = [kron(eye(F),Pm) zeros(dim2(3)*F,(dim1(1)+dim1(2)-2)*F) eye(dim2(3)*F) zeros(dim2(3)*F,(dim2(1)+dim2(2)-2)*F);
    zeros((dim1(1)+dim1(2)-2)*F,dim1(3)*F) eye((dim1(1)+dim1(2)-2)*F) zeros((dim1(1)+dim1(2)-2)*F,dim2(3)*F) -P_tilde];

%CCRB
for s=1:length(sigma_n1)
  CCRB(:,:,s) = ccrb1(CRB1(:,:,s),CRB2(:,:,s),G);
  CCRB_C1(s) = sum(diag(CCRB(1:dim1(3)*F,1:dim1(3)*F,s)));
  CCRB_AB1(s) = sum(diag(CCRB(dim1(3)*F+1:(sum(dim1)-2)*F,dim1(3)*F+1:(sum(dim1)-2)*F,s)));
  CCRB_C2(s) = sum(diag(CCRB((sum(dim1)-2)*F+1:(dim2(3)+sum(dim1)-2)*F,(sum(dim1)-2)*F+1:(dim2(3)+sum(dim1)-2)*F,s)));
  CCRB_AB2(s) = sum(diag(CCRB((dim2(3)+sum(dim1)-2)*F+1:end,(dim2(3)+sum(dim1)-2)*F+1:end,s)));
end

%% Simulations

 Nreal = 200; Ninit = 5; Niter_c = 500; Niter = 500;
P = Nreal*Ninit*length(sigma_n1); p=0; %Progress indicator 

for s=1:length(sigma_n1)
    for n=1:Nreal
        
        %Noisy tensors
        Y1 = X1+sigma_n1(s)*randn(dim1);
        Y2 = X2+sigma_n2(s)*randn(dim2);
        %Initial objective values
        obj = 10^50; obj1 = 10^50; obj2 = 10^50;
        
        for i=1:Ninit
            
            %Initialization
            A1_0 = randn(size(A1)); A1_0(1,:) = ones(1,F);
            B1_0 = randn(size(B1)); B1_0(1,:) = ones(1,F);
            A2_0 = randn(size(A2)); A2_0(1,:) = ones(1,F);
            B2_0 = randn(size(B2)); B2_0(1,:) = ones(1,F);
            C1_0 = randn(size(C1)); C2_0 = randn(size(C2));
            
            % ALS algorithm, uncoupled
            [A1_nu,B1_nu,C1_nu,cost1] = als(Y1,Niter,A1_0,B1_0,C1_0);
            [A2_nu,B2_nu,C2_nu,cost2] = als(Y2,Niter,A2_0,B2_0,C2_0);
            [A2_nu,B2_nu,C2_nu] = amb_correct_1f(A2_nu,B2_nu,C2_nu,Pm*C1_nu);
            
            %Pick best init
            if cost1(end)<obj1
                obj1 = cost1(end);
                A1_u = A1_nu; B1_u = B1_nu; C1_u = C1_nu;
            end
            if cost2(end)<obj2
                obj2 = cost2(end);
                A2_u = A2_nu; B2_u = B2_nu; C2_u = C2_nu;
            end              
            
            % ALS algorithm, coupled
            [A2_n,B2_n,C1_n,cost] = STEREO( Y1,Y2,P1,P2,Pm,Niter_c,(sigma_n2(s)/sigma_n1(s))^2,A2_u,B2_u,C1_u,C2_u, [sigma_n1(s) sigma_n2(s)]);

            %Pick best init
            if cost(end)<obj
                obj = cost(end);
                C1_c = C1_n; A2_c = A2_n; B2_c = B2_n;
                A1_c = P1*A2_c; B1_c = P2*B2_c; C2_c = Pm*C1_c;
            end
            
            p=p+1; clc
            fprintf('Progress %g %%',(p/P)*100)
        end
        
        %Correcting ambiguity
        [A1_u,B1_u,C1_u] = amb_correct(A1_u,B1_u,C1_u,A1,B1,C1);
        [A2_u,B2_u,C2_u] = amb_correct(A2_u,B2_u,C2_u,A2,B2,C2);
        [~,~,C1_c] = amb_correct(A1_c,B1_c,C1_c,A1,B1,C1);
        [A2_c,B2_c,~] = amb_correct(A2_c,B2_c,C2_c,A2,B2,C2);
         
         
        %Manually permute A1_c, B1_c
        mat = A1_c'*A1; ind = [];
        tmp = zeros(size(A1)); tmp2 = zeros(size(B1));
        for i=1:F
              ind(i) = find(mat(:,i)==max(mat(:,i)));
              tmp(:,i) = A1_c(:,ind(i)); tmp2(:,i) = B1_c(:,ind(i));
        end
        A1_c = tmp; B1_c = tmp2;
        
        %Manually permute C2_c
        mat = C2_c'*C2; ind = [];tmp = zeros(size(C2));
        for i=1:F
              ind(i) = find(mat(:,i)==max(mat(:,i)));
              tmp(:,i) = C2_c(:,ind(i)); 
        end
        C2_c = tmp; 
        
       
         
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

end

%% Figures
     
figure
subplot(2,2,1)
semilogy(SNR1,CRB_C1,'k-',SNR1,mse_C1,'r+'); hold on
semilogy(SNR1,CCRB_C1,'k-.',SNR1,mse_C1_c,'bo'); hold on
%semilogy(SNR1,CCRB2_C1,'r-d')
title('Performance bounds for $\mathbf{\theta}_1$','Interpreter','latex')
xlabel('SNR on $\mathcal{X}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
%legend('CRB','CCRB','Non-standard CCRB','Interpreter','latex')
xlim([5 40])

subplot(2,2,2)
semilogy(SNR1,CRB_AB1,'k-',SNR1,mse_AB1,'r+'); hold on
semilogy(SNR1,CCRB_AB1,'k-.',SNR1,mse_AB1_c,'bo'); hold on
%semilogy(SNR1,CCRB2_AB1,'r-d')
title('Performance bounds for $\widetilde{\mathbf{\phi}}_1$','Interpreter','latex')
xlabel('SNR on $\mathcal{X}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
xlim([5 40])

subplot(2,2,3)
semilogy(SNR1,CRB_C2,'k-',SNR1,mse_C2,'r+'); hold on
semilogy(SNR1,CCRB_C2,'k-.',SNR1,mse_C2_c,'bo'); hold on
%semilogy(SNR1,CCRB2_C2,'r-d')
title('Performance bounds for $\mathbf{\theta}_2$','Interpreter','latex')
xlabel('SNR on $\mathcal{X}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
xlim([5 40])

subplot(2,2,4)
semilogy(SNR1,CRB_AB2,'k-',SNR1,mse_AB2,'r+'); hold on
semilogy(SNR1,CCRB_AB2,'k-.',SNR1,mse_AB2_c,'bo'); hold on
%semilogy(SNR1,CCRB2_AB2,'r-d')
title('Performance bounds for $\widetilde{\mathbf{\phi}}_2$','Interpreter','latex')
xlabel('SNR on $\mathcal{X}_1$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
xlim([5 40])

saveas(gcf,'figures/performance_camsap_paper.fig')







