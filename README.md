# Software for Constrained Cramér-Rao bounds with application to hyperspectral super-resolution

Copyright (c) 2019 Clémence Prévost, Konstantin Usevich, Martin Haardt, Pierre Comon, David Brie <br>
Contact: ```clemence.prevost@univ-lorraine.fr```

This MATLAB software reproduces the results from the following article:
```
Incoming
```

## Acknowledgements

Coupled ALS algorithms STEREO and Blind-STEREO are courtesy of C. Kanatsoulis. If you want to use this software for publication, please cite the following paper:
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

## Minimal requirements

- Download and install Tensorlab 3.0: https://www.tensorlab.net

## Demo file
 
 A demo with minimal requirements is available. To proceed, please run the ```demo.m``` file.
 
### Generate the model
 
 ![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%5Cbegin%7Bcases%7D%20%7B%5Cmathcal%7BY%7D%7D_1%20%26%3D%20%5B%5C%21%5B%5Cmathbf%7BA%7D_1%2C%5Cmathbf%7BB%7D_1%2C%5Cmathbf%7BC%7D_1%5D%5C%21%5D%20&plus;%20%5Cmathcal%7BE%7D_1%2C%5C%5C%20%7B%5Cmathcal%7BY%7D%7D_2%20%26%3D%20%5B%5C%21%5B%5Cmathbf%7BA%7D_2%2C%5Cmathbf%7BB%7D_2%2C%5Cmathbf%7BC%7D_2%5D%5C%21%5D%20&plus;%20%5Cmathcal%7BE%7D_2%2C%5C%5C%20%5Cend%7Bcases%7D%20%5Cend%7Balign*%7D)
