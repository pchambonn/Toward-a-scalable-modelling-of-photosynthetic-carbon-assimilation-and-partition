This repository contains the code associated with the publication "Toward a scalable modelling of photosynthetic carbon assimilation and partition during nitrogen starvation for the microalga Chlorella vulgaris NIES 227" by Chambonniere et al. (DOI TO BE ADDED WHEN AVAILABLE).
This code reproduces all the analysis described in the manuscript as well as in the supplementary information provided.

The repository contains the following documents :

  * In the root folder:
    - PhoCAssPaNS_generalCode.RMD
        This document present all the central analysis of the publication and enables to draw the figures shown in the results.
    - Exponential_PI_curve_experiment_massicFit.RMD
        This document aims at fitting the PI curve of the microalgae Chlorella vulgaris NIES 227 at 25°C in exponential growth phase from experiments carried out using the O2 cell set-up.
        Having evidenced the impact of protein quota on the photosynthetic productivity of the microalgae as well as the decrease of cell protein quota under N deficient conditions and under light limited conditions (despite N sufficient conditions) this experiment was designed to test biomass sampled from growth conditions non nutrient starved.
        Measurements were performed at varying light intensities in order to obtain a PI curve and fit the Monod formual modified from Bechet et al. (2015) as presented in the manuscript by Chambonniere et al (2025) to evaluate the growth kinetic parameters Pm and K.
    - PhoCAssPaNS_O2Cell_dataAnalysis.RMD
        This file aims at determining the oxygen productivity of the microalgae CHlorella vulgaris NIES 227 according to light intensity measured as part of the study Scalable modelling of photosynthetic carbon assimilation and partition during
        nitrogen starvation for the microalga Chlorella vulgaris NIES 227. 
        The data collected in that study is loaded and then analysed as described in the manuscript and suplementary information of the associated study.
  
  * In the folder "Auxiliary R functions"
    - o2cell_light_local_interpol.R:
        This function aims at delivering a 3D matrix of light distribution in the Qubit O2 cell set up according to the optical density of the algae solution tested and sofware controlled inicident light intensity.
    - o2cell_PO2_theo_calculation_Pm_massic_lightInhibition_modelModif.R:
        This script aims at computing the O2 productivity in the O2 cell on single experiment knowing the kinetic parameters based on Equation 1 of the manuscript.
    - o2cell_prod_fun_kl_corrected.R:
        The function o2_cell_fun_Kl_corrected is defined to compute algae productivity and respiration rate as measured with the Qubit O2 cell set-up. This calculation takes into account the reaeration rate measured in the O2 cell.
    - pc_interp2.R:
         This function aims at accomplishing bidimensional interpolation.


In order to run, the present folder can be fully downloaded, and the following documents must be added as follows :
  - 
  -
  -
  
