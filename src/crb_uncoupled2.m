function [CRB FIM] = crb_uncoupled2(A1,B1,C1,sigma_n1)

% CRB_UNCOUPLED2 returns uncoupled CRB and FIM for a tensor admitting a 
% CP decomposition[A1,B1,C1] with noise variance sigma_n1
% [CRB FIM] = CRB_UNCOUPLED2(A1,B1,C1,sigma_n1) computes the CRB and FIM.
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

dim1 = [size(A1,1) size(B1,1) size(C1,1)]; F = size(A1,2);

[~,J_13,J_23] = vec_unfold(dim1);
P = [eye(dim1(1)*F) zeros(dim1(1)*F,dim1(2)*F);
     zeros(dim1(2)*F,dim1(1)*F) eye(dim1(2)*F)];
S_C = kron(kr(B1,A1),eye(dim1(3)));
S_A = J_13*kron(kr(C1,B1),eye(dim1(1)));
S_B = J_23*kron(kr(C1,A1),eye(dim1(2)));
D_C = (1/sigma_n1^2)*S_C'*S_C;
D_CAB = (1/sigma_n1^2)*(S_C'*[S_A S_B]*P);
D_AB = (1/sigma_n1^2)*(P'*[S_A'*S_A S_A'*S_B;
                              S_B'*S_A S_B'*S_B]*P);

FIM = [D_C D_CAB;
       D_CAB' D_AB];
                          
CRB = pinv(FIM);
   
end

