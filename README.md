# avo_inversion
Adjoint-state based AVO Inversion Method
Coding Pakage for Computer & Geoscience paper entitled "Constrained Non-linear AVO Inversion based on the Adjoint-State Optimization".
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
