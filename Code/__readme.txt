This readme file is for information and instructions about
how to regenerate all the results in the paper

Table of contents
=======================================================
Regenerate Figure 2 
Regenerate Figure 3
Regenerate Figure 4
Regenerate Figure 5
Regenerate Figure 7
Regenerate Figure 8
Regenerate Figure 9
Regenerate Figure 10
All the filenames and short desciption under folder
=======================================================


===============================================================
#
#  Regenerate Figure 2
#
===============================================================
The script is in "demo_bmc_intro.m"

===============================================================
#
#  Regenerate Figure 3
#
===============================================================
1. Run "Figure3_experiments.m". You only need to change the 
   "FileName" in the last section to the path you want 
   to save the results
2. Run "Figure3_experiments_save.m"
3. Run "Figure3_plot.R"


===============================================================
#
#  Regenerate Figure 4
#
===============================================================
1. Run "Figure4_experiments.m". You only need to change the 
   "FileName" in the last section to the path you want 
   to save the results
2. Run "Figure4_experiments_save.m"
3. Run "Figure4_plot.R" 


===============================================================
#
#  Regenerate Figure 5
#
===============================================================
1. Run "Figure5_experiments.m". You only need to change the 
   "FileName" in the last section to the path you want 
   to save the results
2. Run "Figure5_experiments_save.m"
3. Run "Figure5_plot.R" 



===============================================================
#
#  Regenerate Figure 7
#
===============================================================
1. Run "Figure7_experiments.m". You only need to change the 
   "FileName" in the last section to the path you want 
   to save the results
2. Run "Figure7_experiments_save.m"
3. Run "Figure7_plot.R"


===============================================================
#
#  Regenerate Figure 8
#
===============================================================
1. Run "Figure8_experiments.m". You only need to change the 
   "FileName" in the last section to the path you want 
   to save the results
2. To make plots, follow the same steps as regenerating 
   Figure 3

===============================================================
#
#  Regenerate Figure 9
#
===============================================================
1. Run "Figure9_experiments.m". You only need to change the 
   "FileName" in the last section to the path you want 
   to save the results
2. To make plots, follow the same steps as regenerating 
   Figure 4

===============================================================
#
#  Regenerate Figure 10
#
===============================================================
1. Run "Figure10_experiments.m". You only need to change the 
   "FileName" in the last section to the path you want 
   to save the results
2. To make plots, follow the same steps as regenerating 
   Figure 5




===============================================================
#
#  All the filenames and short desciption under folder
#
===============================================================
*** BMC_CV.m 

A FUNCTION that performs biclustered matrix completion 
with tuning parameters selected by K-fold cross validation on 
MSE

*** BMC_GRID.m

A FUNCTION that performs biclustered matrix completion 
with tuning parameters selected by minmizing BIC using 
grid search

*** BMC_KERNEL_CHECK.m

A FUNCTION to check if S is invertable given X_miss, row and 
column Laplacian

*** BMC_KERNEL_CHECK_TEST.m

A SCRIPT to test BMC_KERNEL_CHECK.m

*** BMC_NM.m

A FUNCTION that performs biclustered matrix completion 
with tuning parameters selected by minmizing BIC using 
nelder-mead algorithm

*** BMC_QN.m

A FUNCTION that solves bmc problem 
using quasi-newton method with log transform

*** BMC_QN_CG.m

A FUNCTION that solves bmc problem using quasi-newton method 
with log transform and use incomplete cholesky and 
conjuate gradient to deal with large and sparse matrix

*** CONFIDENCE_INTERVAL.m

A FUNCTION that computes the confidence interval for the BIC 
exactly

*** DEMO_BMC.m

A SCRIPT showing how to solve bmc problem using different 
methods such as cross-validation, grid search, quasi-newton 
and conjugate-gradient, either with Hutchinson approximation
or not

*** DEMO_BMC_INTRO.m

A SCRIPT to comare matrix completion via singular value 
thresholding and biclustered matrix completion. The comparison
metric is MSE.

*** DEMO_MONDRIAN_EX2.m

A SCRIPT for numerical experiment on the Mondrian data
to compare Quasi-Newton with differnt Hutchinson sample size
and different missing fractions
Comparison metrics inlucdes:
a. BIC 
b. Runtime
c. MSE over missing entries, observed entries and all entries 

*** DEMO_MONDRIAN_EX2_AIC.m

A SCRIPT for numerical experiment on the Mondrian data
to compare Quasi-Newton with differnt Hutchinson sample size
and different missing fractions
Comparison metrics inlucdes:
a. AIC 
b. Runtime
c. MSE over missing entries, observed entries and all entries 

*** Figure3_experiments.m

A SCRIPT to compare between Quasi-Newton with Hutchinson 
estimation (Sample Size = N)
and Quasi-Newton with exact computation (Exact)under 
different missing fractions (0.1, 0.3, 0.5)
Comparison metrics inlucdes:
a. BIC 
b. Runtime
c. MSE over missing entries, observed entries and all entries 

*** Figure3_experiments_save.m

Save the results from Figure3_experiments.m to csv format

*** Figure3_plot.R

Plot the results of Figure3_experiments_save.m using ggplot2

*** Figure4_experiments.m

A SCRIPT to compare among BIC grid search (Grid-N), 
Quasi-Newton with exact computation (Exact) and 
Quasi-Newton with Hutchinson estimation (Sample Size = 5)
under different missing fractions (0.1, 0.3, 0.5)
Comparison metrics inlucdes:
a. BIC 
b. Runtime
c. MSE over missing entries, observed entries and all entries 

*** Figure4_experiments_save.m

Save the results from Figure4_experiments.m to csv format

*** Figure4_plot.R

Plot the results of Figure4_experiments_save.m using ggplot2

*** Figure5_experiments.m

A SCRIPT to compare among cross-validation grid search (Grid-N), 
Quasi-Newton with exact computation (Exact) and Quasi-Newton with 
Hutchinson estimation (Sample Size = 5) under
different missing fractions (0.1, 0.3, 0.5)
Comparison metrics inlucdes:
a. BIC 
b. Runtime
c. MSE over missing entries, observed entries and all entries 

*** Figure5_experiments_save.m

Save the results from Figure5_experiments.m to csv format

*** Figure5_plot.R

Plot the results of Figure3_experiments_save.m using ggplot2

*** Figure7_experiments.m

A SCRIPT for numerical experiment on the Big Radgenomics data
to compare cg with Hutchinson vs Quasi-Neweton with Hutchinson
Comparison metrics inlucdes:
a. BIC 
b. Runtime
c. MSE over missing entries, observed entries and all entries 

*** Figure7_experiments_save.m

Save the results from Figure7_experiments.m to csv format

*** Figure7_plot.R

Plot the results of Figure7_experiments_save.m using ggplot2

*** Figure8_experiments.m

Compared to Figure3_experiments.m, use AIC as objective value

*** Figure9_experiments.m

Compared to Figure4_experiments.m, use AIC as objective value

*** Figure10_experiments.m

Compared to Figure5_experiments.m, use AIC as objective value


*** FMINUNC_PATH.m

A FUNCTION to save the values at each iteration when using 
fminunc to plot the search path

*** GRADIENT_EXACT.m

Compute the gradient of the BIC exactly

*** HESSIAN_EXACT.m

Compute the Hessian of the BIC exactly

*** LAPLACIAN.m

A FUNCTION to generate Laplaican matrix for the simulated data

*** MATRIX_IMPUTE_NESTEROV.m

Impute a set of matrices using Nesterov method

*** MISSING_PATTERN.m

A FUNCTION to impose some missing pattern on a complete matrix 
Patterns include "random", "entire row" or "entire column"

*** NETWORKCOMPONENTS.m

A FUNCTION to find the connected components given an adjacency
matrix. This function is written by Daniel Larremore, 
April 24, 2014, larremor@hsph.harvard.edu

*** SVT.m

Singular value thresholding 

*** TEST_HESSIAN.m

A SCRIPT to test the Hessian and confidence interval computation

*** TUNE_AND_IMPUTE_NESTEROV.m

Impute a set of matrices using Nesterov method