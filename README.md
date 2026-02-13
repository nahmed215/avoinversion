# AVO_Inversion
**Constrained non-linear AVO inversion based on the adjoint-state optimization**

Coding Package for Computer & Geoscience paper entitled "Constrained Non-linear AVO Inversion based on the Adjoint-State Optimization". The pdf copy of this paper is also uploaded (AHMED et al 2022) here, alongwith codes files.

If you use the seismic inversion code, acknowledge it and also cite this article as:

**N. Ahmed, W.W. Weibull, and D. Grana, Constrained non-linear AVO
inversion based on the adjoint-state optimization. Computers and Geosciences (2022), Vol. 168, doi:
https://doi.org/10.1016/j.cageo.2022.105214.**

For the LaTeX editor: 

@article{ahmed2022constrained,\
  title={Constrained non-linear AVO inversion based on the adjoint-state optimization},\
  author={Ahmed, Nisar and Weibull, Wiktor Waldemar and Grana, Dario},\
  journal={Computers \& Geosciences},
  volume={168},
  pages={105214},
  year={2022},
  publisher={Elsevier}
}

This folder contains two main files titled 'Make_Synthetics.m & Non_linear_AVO_Inversion' alongwith exampledata.mat. Keep all three files in a single folder.
To run non-linear avo inversion, the following functions are required (keep all the below functions in another folder with the name Function).
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
First, run 'Make_Synthetics.m', which will create the data sets required for inversion, and then run 'Non_linear_AVO_Inversion'.
