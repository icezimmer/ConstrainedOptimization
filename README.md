## Constrained Optimization using Frank-Wolfe Type Algorithm

We provide the code in MATLAB through the package.
The configurable script 'test.m' is useful to perform and evaluate a Frank-wolfe type algorithm ("Standard" or "Away-Step").
The configurable script 'benchmark.m' compares all the variants of the Frank-Wolfe algorithm selected with all the off-the-shelf algorithm ("interior-point-convex" or "active-set" from quadprog) selected.
The configurable script 'expe_test.m' performs multiple time the Frank-Wolfe type algorithm selected varying the values of the parameters given from the respective lists.
The configurable script 'expe_benchmark.m' compares multiple time all the Frank-Wolfe type algorithms selected and all the off-the-shelf methods selected varying the values of the parameters given from the respective lists.
All the functions codes are in the folder 'src', which we load into the executable scripts via
the command 'add path src'.
The scripts automatically save the results obtained from each experiment in a sub-folder of the folder
'results', named with the date of the experiment.

To test the code and run the experiments, launch 'MATLAB', enter in the main folder,
configure the parameters in one of the configurable scripts and run it.