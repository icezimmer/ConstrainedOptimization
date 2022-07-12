## Constrained Optimization using Frank-Wolfe Method

We provide the code in MATLAB through the package.
The configurable script 'benchmark.m' compares all the variants of the Frank-Wolfe algorithm
and the two off-the-shelf methods provided by the quadratic programming of 'MATLAB'.
Instead, the configurable script 'test.m' is useful to perform and evaluate a precise algorithm.
All the functions codes are in the folder 'src', which we load into the executable scripts via
the command 'add path src'.
The scripts automatically save the results obtained from each experiment in a sub-folder of the folder
'results', named with the date of the experiment.

To test the code and run the experiments, launch 'MATLAB', enter in the main folder,
configure the parameters in one of the configurable scripts and run it in the 'MATLAB' prompt.