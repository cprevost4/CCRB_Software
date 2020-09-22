function [CRB,FIM] = crb_uncoupled(A1,B1,C1,sigma_n1)

% CRB_UNCOUPLED returns uncoupled CRB and FIM for a tensor admitting a 
% CP decomposition[A1,B1,C1] with noise variance sigma_n1
% [CRB FIM] = CRB_UNCOUPLED(A1,B1,C1,sigma_n1) computes the CRB and FIM.
% 
% INPUT ARGUMENTS:
%     A1,B1,C1: low-rank CP factors with rank F
%     sigma_n1: variance of the (white gaussian) noise on the tensor
% OUTPUT ARGUMENTS:
%     CRB: CRB matrix
%     FIM: Fisher information matrix

% Copyright (c) 2020 Clemence Prevost, Konstantin Usevich, Martin Haardt, Pierre Comon, David Brie
% https://github.com/cprevost4/CCRB_Software
% Contact: clemence.prevost@univ-lorraine.fr

%Get dimensions
dim1 = [size(A1,1) size(B1,1) size(C1,1)]; F = size(A1,2);

%Blocks of the FIM
[~,J_13,J_23] = vec_unfold(dim1);
M = eye((dim1(1)+dim1(2))*F);
ind = [1:dim1(1):dim1(1)*F dim1(1)*F+1:dim1(2):(dim1(1)+dim1(2))*F];
M(ind,:) = [];
P = [eye(dim1(1)*F) zeros(dim1(1)*F,dim1(2)*F);
     zeros(dim1(2)*F,dim1(1)*F) eye(dim1(2)*F)];
S_C = kron(kr(B1,A1),eye(dim1(3)));
S_A = J_13*kron(kr(C1,B1),eye(dim1(1)));
S_B = J_23*kron(kr(C1,A1),eye(dim1(2)));
D_C = (1/sigma_n1^2)*S_C'*S_C;
D_CAB = (1/sigma_n1^2)*(S_C'*[S_A S_B]*P*M');
D_AB = (1/sigma_n1^2)*(M*P'*[S_A'*S_A S_A'*S_B;
                              S_B'*S_A S_B'*S_B]*P*M');

%FIM matrix                          
FIM = [D_C D_CAB;
       D_CAB' D_AB];
                          
%Block inversion of the FIM
CRB_AB = inv(D_AB - D_CAB'*inv(D_C)*D_CAB);
CRB_C = inv(D_C) + inv(D_C)*D_CAB*CRB_AB*D_CAB'*inv(D_C);                        
CRB_CAB = -inv(D_C)*D_CAB*CRB_AB;
CRB_ABC = -CRB_AB*D_CAB'*inv(D_C);

%CRB matrix
CRB = [CRB_C CRB_CAB;
       CRB_ABC CRB_AB];

   
end

