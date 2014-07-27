main: Simulations to generate the figures in the report
admm: sequential implementation of the ADMM algorithm in the report
paradmm: parallel implemetation using Matlab parfor. This is very basic implementation with suboptimal performance.

To generate the results:
Run the main file as is to generate the figures in the report. The communication overhead is 
computed using the given formulas in the report.

To run the parallel implementation:
1. Uncomment line 43
2. comment line 46
3. Properly configure matlab parallel local profile J=number of core