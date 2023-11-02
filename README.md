# Constrained Optimization using Frank-Wolfe Type Algorithm

This repository provides MATLAB code for performing constrained optimization using Frank-Wolfe type algorithms. The code is organized into a package, and several configurable scripts are provided to facilitate different types of analyses and comparisons.

## Available Scripts

### `test.m`
This script allows you to perform and evaluate a Frank-Wolfe type algorithm, either "Standard" or "Away-Step". You can configure the script to suit your specific needs.

### `benchmark.m`
This script enables comparison between all selected variants of the Frank-Wolfe algorithm and selected off-the-shelf algorithms from `quadprog` ("interior-point-convex" or "active-set"). 

### `expe_test.m`
This script performs multiple runs of a selected Frank-Wolfe type algorithm, varying the parameter values based on provided lists. 

### `expe_benchmark.m`
This script conducts multiple comparisons between all selected Frank-Wolfe type algorithms and off-the-shelf methods, varying the parameter values based on provided lists.

## Code Structure

All function codes are located in the `src` folder. To use these functions in the executable scripts, add the source folder to the MATLAB path using the command `addpath src`. 

Results from each experiment are automatically saved in a sub-folder within the `results` folder. The sub-folder is named based on the date of the experiment.

## Getting Started

To run the experiments and test the code:

1. Launch MATLAB and navigate to the main folder of this repository.
2. Open one of the configurable scripts (`test.m`, `benchmark.m`, `expe_test.m`, or `expe_benchmark.m`).
3. Configure the parameters within the script as needed.
4. Run the script.

Explore and modify the scripts as needed to suit your specific use case. Enjoy optimizing!
