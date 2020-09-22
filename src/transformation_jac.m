function [H2, H_CCRB] = transformation_jac(A2,B2,C1,P1,P2,Pm,alpha,beta,dim1,dim2)

% TRANSFORMATION_JAC returns matrices containing the Jacobian of the
% non-linear constraints, according to their definition in eq.(20).
% [H2, H_CCRB] = TRANSFORMATION_JAC(A2,B2,C1,P1,P2,Pm,alpha,beta,dim1,dim2) computes the Jacobian.
% 
% INPUT ARGUMENTS:
%     A2,B2,C1: low-rank CP factors with rank F
%     P1,P2,Pm: degradation matrices
%     alpha,beta: scaling factors for non-linear constraints
%     dim1,dim2: dimensions of the two low-resolution images
% OUTPUT ARGUMENTS:
%     H2,H_CCRB: Jacobian matrices. They contain the same elements but
%     ordered differently. For the "classic" CCRB, one should use H_CCRB.

% Copyright (c) 2020 Clemence Prevost, Konstantin Usevich, Martin Haardt, Pierre Comon, David Brie
% https://github.com/cprevost4/CCRB_Software
% Contact: clemence.prevost@univ-lorraine.fr

%Get dimensions
F = size(A2,2);

%Partial derivatives
dA1_dA2 = kron(diag(1./alpha),P1(2:end,2:end)) - ...
    (kron(diag(1./(alpha.^2)),P1(2:end,2:end))*kr(eye(F),A2(2:end,:))*kron(eye(F),P1(1,2:end)));
dB1_dB2 = kron(diag(1./beta),P2(2:end,2:end)) - ...
    (kron(diag(1./(beta.^2)),P2(2:end,2:end))*kr(eye(F),B2(2:end,:))*kron(eye(F),P2(1,2:end)));
dC2_dC1 = kron(diag(1./(alpha.*beta)),Pm);
dC2_dA2 = -kron(diag(1./(beta.*(alpha.^2))),Pm)*kr(eye(F),C1)*kron(eye(F),P1(1,2:end));
dC2_dB2 = -kron(diag(1./(alpha.*(beta.^2))),Pm)*kr(eye(F),C1)*kron(eye(F),P2(1,2:end));

%Jacobian for reparametrized CRB
H2 = blkdiag(dA1_dA2,dB1_dB2,dC2_dC1);
H2((dim1(1)+dim1(2)-2)*F+1:end,1:(dim2(1)-1)*F) = dC2_dA2;
H2((dim1(1)+dim1(2)-2)*F+1:end,(dim2(1)-1)*F+1:(dim2(1)+dim2(2)-2)*F) = dC2_dB2;

%Jacobian for standard CCRB
H_CCRB = [-dC2_dC1 zeros(dim2(3)*F,(dim1(1)+dim1(2)-2)*F) eye(dim2(3)*F) -dC2_dA2 -dC2_dB2;
          zeros((dim1(1)+dim1(2)-2)*F,dim1(3)*F) eye((dim1(1)+dim1(2)-2)*F)...
                            zeros((dim1(1)+dim1(2)-2)*F,dim2(3)*F) -blkdiag(dA1_dA2,dB1_dB2)];

end

