# Software for Constrained Cramér-Rao bounds with application to hyperspectral super-resolution

Copyright (c) 2020 Clémence Prévost, Konstantin Usevich, Martin Haardt, Pierre Comon, David Brie <br>
Contact: ```clemence.prevost@univ-lorraine.fr```

This MATLAB software reproduces the results from the following papers:

```
@inproceedings{prevost2019performance,
  title={Performance bounds for coupled CP model in the framework of hyperspectral super-resolution},
  author={Prevost, Cl{\'e}mence and Usevich, Konstantin and Comon, Pierre and Haardt, Martin and Brie, David},
  booktitle={2019 IEEE 8th International Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP)},
  pages={201--205},
  year={2019},
  organization={IEEE}
}
```

```
Incoming
```

## Acknowledgements

Uncoupled ALS algorithm is courtesy of R. Bro.
Coupled ALS algorithms STEREO and Blind-STEREO are courtesy of C. Kanatsoulis. If you want to use this software for publication, please cite the following papers:
```
@inproceedings{bro1998multi,
  title={Multi-way analysis in the food industry-models, algorithms, and applications},
  author={Bro, Rasmus},
  booktitle={MRI, EPG and EMA,” Proc ICSLP 2000},
  year={1998},
  organization={Citeseer}
}
```

```
@article{kanatsoulis2018hyperspectral,
  title={Hyperspectral super-resolution: A coupled tensor factorization approach},
  author={Kanatsoulis, Charilaos I and Fu, Xiao and Sidiropoulos, Nicholas D and Ma, Wing-Kin},
  journal={IEEE Transactions on Signal Processing},
  volume={66},
  number={24},
  pages={6503--6517},
  year={2018},
  publisher={IEEE}
}
``` 

## Content

 - demo.m : demo file with minimal requirements 
 -/data : contains simulations results, for speed purpose. These results can also be obtained in the demo file.
 - /demos : contains demo files that produce tables and figures (including ```main.m```

 - /figures : where the figures are saved
 - /images : contains illustrative figures for this ```README.md```
 - /src : contains helpful files to run the demos

## Minimal requirements

In order to run the demo file and reproduce the figures, you will need to:
- Download and install Tensorlab 3.0: https://www.tensorlab.net

## Demo file
 
 A demo with minimal requirements is available. To proceed, please run the ```demo.m``` file.
 
### Generate the model
 
 ![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%5Cbegin%7Bcases%7D%20%7B%5Cmathcal%7BY%7D%7D_1%20%26%3D%20%5B%5C%21%5B%5Cmathbf%7BA%7D_1%2C%5Cmathbf%7BB%7D_1%2C%5Cmathbf%7BC%7D_1%5D%5C%21%5D%20&plus;%20%5Cmathcal%7BE%7D_1%2C%5C%5C%20%7B%5Cmathcal%7BY%7D%7D_2%20%26%3D%20%5B%5C%21%5B%5Cmathbf%7BA%7D_2%2C%5Cmathbf%7BB%7D_2%2C%5Cmathbf%7BC%7D_2%5D%5C%21%5D%20&plus;%20%5Cmathcal%7BE%7D_2%2C%5C%5C%20%5Cend%7Bcases%7D%20%5Cend%7Balign*%7D%20%5Cbegin%7Balign*%7D%20%5Cquad%5Ctext%7B%20s.%20to%20%7D%5Cmathbf%7BA%7D_1%20%3D%20%5Cmathbf%7BP%7D%5Cmathbf%7BA%7D_2%5Ccdot%5Cmathbf%7B%5Calpha%7D%5E%7B-1%7D%2C%20%5Cmathbf%7BB%7D_1%20%3D%20%5Cmathbf%7BQ%7D%5Cmathbf%7BB%7D_2%5Ccdot%5Cmathbf%7B%5Cbeta%7D%5E%7B-1%7D%2C%5Cnonumber%5C%5C%20%5Cmathbf%7BC%7D_2%20%3D%5Cmathbf%7BR%7D%5Cmathbf%7BC%7D_1%5Ccdot%28%5Cmathbf%7B%5Calpha%5Cbeta%7D%29%5E%7B-1%7D%5Cnonumber%20%5Cend%7Balign*%7D)
 
 ### Calculate the Cramér-Rao bounds
 
 We compute the following lower bounds:
 - Uncoupled CRB;
 - Fully-coupled CCRB;
 - Blind-CCRB;
 
 for the following parameters:
 
 ![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%7B%5Cwidetilde%7B%5Cmathbf%7B%5Ctheta%7D%7D%7D%5E%7BT%7D%20%26%3D%20%5Cbegin%7Bbmatrix%7D%5Ctext%7Bvec%7D%28%7B%28%5Cmathbf%7BA%7D_2%29_%7B2%3AI%2C%3A%7D%7D%29%5E%7BT%7D%20%26%20%5Ctext%7Bvec%7D%28%7B%28%5Cmathbf%7BB%7D_2%29_%7B2%3AJ%2C%3A%7D%7D%29%5E%7BT%7D%20%26%20%5Ctext%7Bvec%7D%28%5Cmathbf%7BC%7D_1%29%5E%7BT%7D%5Cend%7Bbmatrix%7D%2C%20%5C%5C%20%7B%5Cwidetilde%7B%5Cmathbf%7B%5Cphi%7D%7D%7D%5E%7BT%7D%20%26%3D%20%5Cbegin%7Bbmatrix%7D%5Ctext%7Bvec%7D%28%7B%28%5Cmathbf%7BA%7D_1%29_%7B2%3AI_H%2C%3A%7D%7D%29%5E%7BT%7D%20%26%20%5Ctext%7Bvec%7D%28%7B%28%5Cmathbf%7BB%7D_1%29_%7B2%3AJ_H%2C%3A%7D%7D%29%5E%7BT%7D%20%26%20%5Ctext%7Bvec%7D%28%5Cmathbf%7BC%7D_2%29%5E%7BT%7D%5Cend%7Bbmatrix%7D.%20%5Cend%7Balign*%7D)
 
  ### Run the algorithms
  
  We run the following algorithms:
  - Uncoupled ALS;
  - STEREO;
  - Blind-STEREO;
  
  on our model.
  The number of realizations, as well as the maximum number of iterations, can be user-specified in the ```Pre-allocation``` section of the file.
  We compute the MSE provided by each algorithm and compare it to the bounds.

  ### Display the results 
  
  We plot the bounds and total MSEs as a function of the noise on the first tensor (HSI).
  
  <img src="images/results.jpg?raw=true"/>
  
  ## Reproduce figures from the paper
  
  To do so, you need to run the ```main.m``` file. Here, a menu is available and allows you to choose which figure or table you want to generate. Each number in the table below corresponds to a set of figures.

