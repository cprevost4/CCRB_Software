function [CRB, D] = crb_allatonce(A1,B1,C1,A2,B2,C2,sigma_n1,sigma_n2)

% CRB_ALLATONCE returns uncoupled CRB and FIM for a tensor admitting a 
% CP decomposition[A1,B1,C1] with noise variance sigma_n1
% [CRB D] = CRB_ALLATONCE(A1,B1,C1,A2,B2,C2,sigma_n1,sigma_n2) computes the CRB and FIM.
% 
% INPUT ARGUMENTS:
%     A1,B1,C1,A2,B2,C2: low-rank CP factors with rank F
%     sigma_n1,sigma_n2: variance of the (white gaussian) noise on the
%     tensors
% OUTPUT ARGUMENTS:
%     CRB: CRB matrix
%     FIM: Fisher information matrix

% Copyright (c) 2020 Clemence Prevost, Konstantin Usevich, Martin Haardt, Pierre Comon, David Brie
% https://github.com/cprevost4/CCRB_Software
% Contact: clemence.prevost@univ-lorraine.fr

%Extract dimensions
dim1 = [size(A1,1) size(B1,1) size(C1,1)]; 
dim2 = [size(A2,1) size(B2,1) size(C2,1)]; F = size(A1,2);
%Permutation matrices
[J_12_1,J_13_1,~] = vec_unfold(dim1);
[J_12_2,J_13_2,~] = vec_unfold(dim2);
%Covariance matrix
Sigma = blkdiag((1/sigma_n1^2)*eye(prod(dim1)),(1/sigma_n2^2)*eye(prod(dim2)));

%Selection matrix
ind = [1:dim2(1):dim2(1)*F dim2(1)*F+1:dim2(2):(dim2(1)+dim2(2))*F];
M1 = eye((dim2(1)+dim2(2)+dim1(3))*F); M1(ind,:) = []; %For omega
ind = [1:dim1(1):dim1(1)*F dim1(1)*F+1:dim1(2):(dim1(1)+dim1(2))*F];
M2 = eye((dim1(1)+dim1(2)+dim2(3))*F); M2(ind,:) = []; %For epsilon
M = blkdiag(M1,M2);

%Tensor unfoldings
S_A1 = kron(kr(C1,B1),eye(dim1(1)));
S_B1 = J_12_1'*kron(kr(C1,A1),eye(dim1(2)));
S_C1 = J_13_1'*kron(kr(B1,A1),eye(dim1(3)));
S_A2 = kron(kr(C2,B2),eye(dim2(1)));
S_B2 = J_12_2'*kron(kr(C2,A2),eye(dim2(2)));
S_C2 = J_13_2'*kron(kr(B2,A2),eye(dim2(3)));

%Derivative
Mu = [zeros(size(S_C1,1),size(S_A2,2)+size(S_B2,2)) S_C1 S_A1 S_B1 zeros(size(S_C1,1),size(S_C2,2));
      S_A2 S_B2 zeros(size(S_A2,1),size(S_C1,2)+size(S_A1,2)+size(S_B1,2)) S_C2];

%FIM and CRB
D = M*Mu'*Sigma*Mu*M';
CRB = pinv(D);
end

