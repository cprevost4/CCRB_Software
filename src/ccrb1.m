function CCRB = ccrb1(CRB1,CRB2,G)

% CCRB1 returns the CCRB as expressed by Gorman and Hero
% CCRB = CCRB1(CRB1,CRB2,G) computes the CCCRB.
% 
% INPUT ARGUMENTS:
%     CRB1,CRB2: CRB matrices obtained from uncoupled tensors 
%     G: Jacobian of the constraints on the parameters
% OUTPUT ARGUMENTS:
%     CCRB: CCRB matrix

% Copyright (c) 2020 Clemence Prevost, Konstantin Usevich, Martin Haardt, Pierre Comon, David Brie
% https://github.com/cprevost4/CCRB_Software
% Contact: clemence.prevost@univ-lorraine.fr

CRB = [CRB1 zeros(size(CRB1,1),size(CRB2,2));
              zeros(size(CRB2,1),size(CRB1,2)) CRB2];
CCRB = CRB - CRB*G'*inv(G*CRB*G')*G*CRB;

end

