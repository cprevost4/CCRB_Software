%                Constrained Cramér-Rao lower bounds for                  %
%         hyperspectral super-resolution: parametrization and             %
%                design of the system of measurements                     %
%-------------------------------------------------------------------------%

% Copyright (c) 2020 Clemence Prevost, Konstantin Usevich, Martin Haardt,
% Pierre Comon, David Brie
% https://github.com/cprevost4/CCRB_Software
% Contact: clemence.prevost@univ-lorraine.fr

% This software reproduces the results from the papers called:
% (1) "Performance bounds for coupled CP model in the framework of 
% hyperspectral super-resolution - C.Prévost, K.Usevich, M. Haardt, 
% D.Brie, P.Comon.
% (2) "Constrained Cramér-Rao lower bounds for hyperspectral super-resolution:
% parametrization and design of the system of measurements - C.Prévost, 
% K.Usevich, M. Haardt, D.Brie, P.Comon.

% In order to run the demo, you will need to add to your MATLAB path:
% - Tensorlab 3.0: https://www.tensorlab.net
% - Codes by C. Kanatsoulis: 
%      https://github.com/marhar19/HSR_via_tensor_decomposition
% - Codes for HySure: https://github.com/alfaiate/HySure
% - Codes for FUSE: https://github.com/qw245/BlindFuse
% - Hyperspectral data :
%      https://github.com/marhar19/HSR_via_tensor_decomposition

%-------------------------------------------------------------------------%
%                              CONTENT
% - /data : contains data for synthetic examples (Section VI.D)
% - /demos : contains demo files that produce tables and figures
% - /figures : where the tables and figures are saved
% - /src : contains helpful files to run the demos
%-------------------------------------------------------------------------%
%                                MENU
% You can launch a specified demo by typing its number. The resulting tables
% and figures produced will be stored in the figures folder.
%
% 1:  produces Figure 1 from (1) 
% 2:  produces Fig. 1 and 3 from (2) 
% 3:  produces Fig. 2 and 4 from (2) 
% 4:  produces Fig. 6 and 7 from (2) 
% 5:  produces Fig. 9 and 10 from (2)
% 6:  produces Fig. 11 and 12 from (2)
% 7:  produces Fig. 16, 18 or 20 from (2) 
% 8:  produces Fig. 15, 17 or 19 from (2) 

%-------------------------------------------------------------------------%

list_demos = ["camsap_paper" "scenario1_performance" "scenario2_performance" ...
              "choice_rank" "influence_q" "influence_d" ...
              "choice_lambda_scenario2" "choice_lambda_scenario1"
    ];

prompt = "Which file do you want to run ?";
num = input(prompt);
eval(list_demos(num));