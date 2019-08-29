## lawic: Large Amplitude water Waves Interacting with underlying Currents

Lawic solves large amplitude steady periodic water waves traveling on rotational flow. At present, it is distrbuted as an executable program with input parameters that can be changed. 

### Get started

Simply find the excutable file in the folder _lawic-1.0.0-macos64/_ for Apple user or in the folder _lawic-1.0.0-win64/_ for windows user and run _lawic_ along with the input file params.dat. Parameters in the input file are self-explanatory. Note that two types of problems are solved respectively by setting:
1. negative relative mass flux p0, and 
2. positive relative mass flux p0. 

For the first case, nonzero vorticity is valid while for the second case the value of the vorticity will be also taken as zero no matter what is provided in the input file. When solving the first type of problem, an output file named "fix_p0.txt" will be generated; when solving the second type of problem, a file named "large_u.txt" will be generated.

The output files can then be loaded and processed using the tools in the folder _postproc/_. There are two files in the folder:
1. flow.py which is a class file with all the functions for postprocessing.
2. fix_p0.py which is an example for using flow.py to load the data and use the functions for generating figures.

### Mathematics

Under consideration is two-dimensional steady periodic traveling surface waves propagating over water of a finite depth. Flat bed is assumed. Both irrotational and rotational flows can be considered. In this version, only rotational flow with constant vorticity is solved. The motion of inviscid fluids under gravity is governed by Euler's equation with kinematic surface and dynamic boundary conditions at the surface and the kinematic boundary condition at the bed which is assumed to be impenetrable. For more on the mathematics and the proof of the existence of large amplitude wave solutions to this problem, refer to [5].

For solving large amplitude waves, the free surface problem is first transformed to a fixed boundary problem by using the Dubreril-Jacotin (DJ)transformation. Secondly, the rectangular domain is discretized and the equation is discretized using the finite difference method. Center difference is used at intermediate grid points and backward/forward finite difference is used at boundary points. The nonlinear algebriac equations are then solved using continuation method from the laminar flow solution to large amplitude waves close to waves with stagnation points. See [2-5] for details of the numerical methods.

Most of the relevant studies focus on the case of _u_<_c_ and hence p0<0, while in [6] it has proved the existence of the water waves with _u_>_c_ (p0>0) for irrotational flow. The latter case is also solved. 


### Examples
The 

### References

1. Chen, L. & Basu, B. (2019). [Numerical investigations of two-dimensional irrotational water waves over finite-depth with uniform current](https://www.tandfonline.com/doi/full/10.1080/00036811.2019.1636974). _Applicable Analysis_. 
1. Ko, J. & Strauss, W. (2008). [Large-amplitude steady rotational water waves](https://doi.org/10.1016/j.euromechflu.2007.04.004). _European Journal of Mechanics / B Fluids_, 27(2), 96-109. 
1. Amann, D. & Kalimeris, K. (2018). [A numerical continuation approach for computing water waves of large wave height](https://doi.org/10.1016/j.euromechflu.2017.10.001). _European Journal of Mechanics / B Fluids_, 67, 314-328. 
1. Constantin, A. & Strauss, W. (2004). [Exact steady periodic water waves with vorticity](https://onlinelibrary.wiley.com/doi/pdf/10.1002/cpa.3046). _Communications on Pure Applied Mathematics_, 57(4), 481-527.
1. Constantin, A. (2011). _Nonlinear water waves with applications to wave-current interactions and tsunamis_, Vol. 81 of CMBS-NSF Reg. Conf. Ser. Appl. Math., SIAM Philadelphia.
1. Basu, B. (2018). [On the existence of two-dimensional irrotational water waves over finite depth with uniform current](https://doi.org/10.1080/00036811.2017.1376321). _Applicable Analysis_, 97(14), 2523-2532.

### License

### Authors
Dr. [Lin Chen](https://chen-lin.github.io) <br/>
Prof. [Biswajit Basu](https://www.tcd.ie/research/profiles/?profile=basub)
