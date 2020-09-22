function [ A_hat,B_hat,C,A_tilde,B_tilde,C_tilde,cost ] = Blind_STEREO(H,M,P3,MAXIT,A_hat,B_hat,A_tilde,B_tilde,C_tilde,lamda,sigma_n)
% Blind STEREO algorithm
% (c) Charilaos I. Kanatsoulis, University of Minnesota, Jan 7 , 2018
% nikos@umn.edu
% 
% Reference 1: C.I. Kanatsoulis, X. Fu, N.D. Sidiropoulos and W.K. Ma, 
%``Hyperspectral Super-resolution: A Coupled Tensor Factorization
%Approach,'' IEEE Transactions in Signal Processing, 2018

% Reference 2: C.I. Kanatsoulis, X. Fu, N.D. Sidiropoulos and W.K. Ma, 
%``Hyperspectral Super-resolution via Coupled Tensor Factorization:
%Identifiability and Algorithms,'' IEEE International Conference on 
%Acoustics, Speech and Signal Processing (ICASSP), 2018

[Ih,Jh,Kh]=size(H); K=Kh;
[I,J,Km]=size(M);
%% create the matrix equivalent models for the input tensor
% H1=zeros(Kh*Jh,Ih);
% H2=zeros(Kh*Ih,Jh);
% H3=zeros(Ih*Jh,Kh);

H1=reshape(H,[Ih,Jh*Kh])';

H2=reshape(permute(H,[2 1 3]),[Jh,Kh*Ih])';

H3=reshape(H,[Ih*Jh,Kh]);

nH=norm(H3,'fro');

M1=reshape(M,[I,J*Km])';

M2=reshape(permute(M,[2 1 3]),[J,Km*I])';

M3=reshape(M,[I*J,Km]);


nM=norm(M3,'fro');
%% initialize
temp_tilde=kr(B_tilde,A_tilde);
C=(temp_tilde\H3)';


%% Alternating least squares
eps=1e-03;
cost(1)=inf; diff_cost = inf; iter = 1;
temp3m=kr(B_hat,A_hat);
Km= (B_hat'*B_hat).*(A_hat'*A_hat);
inv_K=pinv(Km);
M_ferror=norm(M3-temp3m*C_tilde','fro');

while (iter < MAXIT) && (diff_cost > eps) && (cost(iter) > eps)
    
    iter = iter+1;
    %     iter
    temp1h=kr(C,B_tilde);
    Kh= lamda*(C'*C).*(B_tilde'*B_tilde);
    A_tilde=(Kh\(lamda*temp1h'*H1))';
    %norm_At = sqrt(sum(A_tilde.^2)); A_tilde = A_tilde.*repmat(1./norm_At,Ih,1);
    
    temp2h=kr(C,A_tilde);
    Kh= lamda*(C'*C).*(A_tilde'*A_tilde);
    B_tilde=(Kh\(lamda*temp2h'*H2))';
    %norm_Bt = sqrt(sum(B_tilde.^2)); B_tilde = B_tilde.*repmat(1./norm_Bt,Jh,1);
    
    
%     temp3h=kr(B_tilde,A_tilde);
%     Kh=lamda*(B_tilde'*B_tilde).*(A_tilde'*A_tilde);
%     As=P3'*P3;
%     Bs=Kh*inv_K;
%     Cs=(lamda*H3'*temp3h+P3'*M3'*temp3m)*inv_K;
%     C=sylvester(full(As),Bs,Cs);
%     C_tilde=P3*C;
%     
%         S1_hat1=khatri_rao(C,B_hat)*A_hat';
%     norm(S1-S1_hat1,'fro')/norm(S1,'fro')
    
    
    temp1m=kr(C_tilde,B_hat);
    Kh= (C_tilde'*C_tilde).*(B_hat'*B_hat);
    A_hat=(Kh\(temp1m'*M1))';
    %norm_A = sqrt(sum(A_hat.^2)); A_hat = A_hat.*repmat(1./norm_A,I,1);
    
    temp2m=kr(C_tilde,A_hat);
    Kh= (C_tilde'*C_tilde).*(A_hat'*A_hat);
    B_hat=(Kh\(temp2m'*M2))';
    %norm_B = sqrt(sum(B_hat.^2)); B_hat = B_hat.*repmat(1./norm_B,J,1);
    
    
    temp3m=kr(B_hat,A_hat);
    Km= (B_hat'*B_hat).*(A_hat'*A_hat);
    inv_K=pinv(Km);
    temp3h=kr(B_tilde,A_tilde);
    Kh=lamda*(B_tilde'*B_tilde).*(A_tilde'*A_tilde);
    As=P3'*P3;
    Bs=Kh*inv_K;
    Cs=(lamda*H3'*temp3h+P3'*M3'*temp3m)*inv_K;
    C=sylvester(full(As),Bs,Cs);
    C_tilde=P3*C;
    %     %     C_hat=P3*C;
    %         C1=(temp3h\H3t)';
%     S1_hat1=khatri_rao(C,B_hat)*A_hat';
%     norm(S1-S1_hat1,'fro')/norm(S1,'fro')
    %     S1_hat1=khatri_rao(V*C,B)*A';
    %     norm(S1-S1_hat1,'fro')/norm(S1,'fro')
        M_ferror=1/(sigma_n(2)^2)*norm(M3-temp3m*C_tilde','fro');
        H_ferror=1/(sigma_n(1)^2)*norm(H3-temp3h*C','fro');
        cost(iter)=lamda*H_ferror+M_ferror;
    diff_cost = cost(iter-1)-cost(iter);
    %     %     drawnow;semilogy(cost);
    %     if (H_ferror/nH)<eps && (M_ferror/nM)<eps
    %         break
    %     end
end

C =     C.*repmat(A_tilde(1,:).*B_tilde(1,:),K,1);
C_tilde =     C_tilde.*repmat(A_hat(1,:).*B_hat(1,:),size(M,3),1);
A_hat =     A_hat.*repmat(1./A_hat(1,:),I,1);
B_hat =     B_hat.*repmat(1./B_hat(1,:),J,1);
A_tilde =     A_tilde.*repmat(1./A_tilde(1,:),Ih,1);
B_tilde =     B_tilde.*repmat(1./B_tilde(1,:),Jh,1);

end
