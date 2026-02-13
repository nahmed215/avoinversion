# avoinversion
*Constrained non-linear AVO inversion based on the adjoint-state optimization*

Coding Pakage for Computer & Geoscience paper entitled "Constrained Non-linear AVO Inversion based on the Adjoint-State Optimization". The pdf copy of this paper is also uploaded (AHMED et al 2022) here alongwith codes files.

Please cite this article as: N. Ahmed, W.W. Weibull and D. Grana, Constrained non-linear AVO
inversion based on the adjoint-state optimization. Computers and Geosciences (2022), Vol. 168, doi:
https://doi.org/10.1016/j.cageo.2022.105214.

This folder contains two main files titled 'Make_Synthetics.m & Non_linear_AVO_Inversion' alongwith exampledata.mat. Kepp all three files in a single folder.
To run non-linear avo inversion, following function are required (keep all the below functions in other folder with name Function).
Functions:
\Functions\Aki_Richards.m                           
\Functions\cubic_min1.m                            
\Functions\cubic_min2.m                         
\Functions\fwmod.m                                  
\Functions\linesearch_backtrack.m                   
\Functions\linesearch_more.m                        
\Functions\quad_min1.m                             
\Functions\quad_min2.m
\Functions\Ricker.m                             
\Functions\redblue.m                               
\Functions\update_trial_int.m
\Functions\vel_smoother.m                      
\Functions\wei_lbfgs_hess.m                       
\Functions\wiggle.m
First, run 'Make_Synthetics.m' it will create a data sets required for inversion and then run 'Non_linear_AVO_Inversion'.
