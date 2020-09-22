function CCRB = ccrb3(A2,B2,C1,dim1,dim2,sigma_n1,sigma_n2,P1,P2,Pm)

% CCRB3 returns the reparametrized CRB 
% CCRB = CCRB3(A2,B2,C1,dim1,dim2,sigma_n1,sigma_n2,P1,P2,Pm) computes the
% reparametrized CRB.
% 
% INPUT ARGUMENTS:
%     A2,B2,C1: CP factors of the SRI
%     dim1,dim2: dimensions of the HSI and MSI
%     sigma_n1,sigma_n2: noise level of the HSI and MSI
%     P1,P2,Pm: degradation matrices
% OUTPUT ARGUMENTS:
%     CCRB: Reparametrized CRB matrix

% Copyright (c) 2020 Clemence Prevost, Konstantin Usevich, Martin Haardt, Pierre Comon, David Brie
% https://github.com/cprevost4/CCRB_Software
% Contact: clemence.prevost@univ-lorraine.fr

F = size(A2,2); 
alpha = P1(1,:)*A2; beta = P2(1,:)*B2;

[J_12_1,J_13_1,~] = vec_unfold([dim1(1) dim1(2) dim1(3)]);
[J_12_2,J_13_2,~] = vec_unfold([dim2(1) dim2(2) dim2(3)]);
ind = [1:dim2(1):dim2(1)*F dim2(1)*F+1:dim2(2):(dim2(1)+dim2(2))*F];
M = eye((dim2(1)+dim2(2)+dim1(3))*F); M(ind,:) = [];

Sigma = blkdiag((1/sigma_n1^2)*eye(prod(dim1)*F),(1/sigma_n2^2)*eye(prod(dim2)*F));

S_A1 = kron(kr(eye(F),kr(C1,P2*B2/diag(beta))),P1)*...
    (kron(diag(1./alpha),eye(dim2(1))) - kron(diag(1./alpha.^2),eye(dim2(1)))*kr(eye(F),A2)*kron(eye(F),P1(1,:)));
S_B1 = kron(eye(F),J_12_1')*kron(kr(eye(F),kr(C1,P1*A2/diag(alpha))),P2)*...
    (kron(diag(1./beta),eye(dim2(2))) - kron(diag(1./beta.^2),eye(dim2(2)))*kr(eye(F),B2)*kron(eye(F),P2(1,:)));
S_C1 = kron(eye(F),J_13_1')*kron(kr(eye(F),kr(P2*B2/diag(beta),P1*A2/diag(alpha))),eye(dim1(3)));

S_A2 = kron(kr(eye(F),kr(Pm*C1/diag(beta),B2)),eye(dim2(1)))*...
    (kron(diag(1./alpha),eye(dim2(1)))- kron(diag(1./alpha.^2),eye(dim2(1)))*kr(eye(F),A2)*kron(eye(F),P1(1,:)));
S_B2 = kron(eye(F),J_12_2')*kron(kr(eye(F),kr(Pm*C1/diag(alpha),A2)),eye(dim2(2)))*...
    (kron(diag(1./beta),eye(dim2(2)))- kron(diag(1./beta.^2),eye(dim2(2)))*kr(eye(F),B2)*kron(eye(F),P2(1,:)));
S_C2 = kron(eye(F),J_13_2')*kron(kr(eye(F),kr(B2,A2)),Pm)*kron(diag(1./(alpha.*beta)),eye(dim1(3)));

Mu = [S_A1 S_B1 S_C1;
      S_A2 S_B2 S_C2];

D = M*Mu'*Sigma*Mu*M';

CCRB = inv(D);


end

