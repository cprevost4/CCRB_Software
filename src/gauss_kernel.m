function Phi = gauss_kernel(q,sigma)

% GAUSS_KERNEL returns Gaussian blurring kernel along one dimension
% Phi = GAUSS_KERNEL(q,sigma) computes kernel Phi of size 1xq with parameter sigma
% 
% INPUT ARGUMENTS:
%     q : size of kernel 
%     sigma: normal distribution parameter (default: sigma = 0.5)
% OUTPUT ARGUMENTS:
%     Phi: Gaussian kernel of size 1xq

% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Software
% Contact: clemence.prevost@univ-lorraine.fr

if nargin==1
    sigma = q/4/sqrt(2*log(2));
end

x = linspace(-q / 2, q / 2, q);
Phi = exp(-x .^ 2 / (2* sigma ^ 2));

end

