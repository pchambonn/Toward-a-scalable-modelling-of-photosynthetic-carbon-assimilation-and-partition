This repository contains the code associated with the publication "_Toward a scalable modelling of photosynthetic carbon assimilation and partition during nitrogen starvation for the microalga Chlorella vulgaris NIES 227_" by Chambonniere et al. (DOI TO BE ADDED WHEN AVAILABLE).
This code reproduces all the analysis described in the manuscript as well as in the supplementary information provided.

The repository contains the following documents :

  * In the root folder:
    - **PhoCAssPaNS_generalCode.RMD** [R markdown document]
        This document present all the central analysis of the publication and enables to draw the figures shown in the results.
    - **Exponential_PI_curve_experiment_massicFit.RMD** [R markdown document]
        This document aims at fitting the PI curve of the microalgae _Chlorella vulgaris_ NIES 227 at 25°C in exponential growth phase from experiments carried out using the O2 cell set-up.
    - PhoCAssPaNS_O2Cell_dataAnalysis.RMD [R markdown document]
        This file aims at determining the oxygen productivity of the microalgae _Chlorella vulgaris_ NIES 227 according to light intensity measured as part of the study '_Scalable modelling of photosynthetic carbon assimilation and partition during nitrogen starvation for the microalga Chlorella vulgaris NIES 227_'. 
        The data collected in that study is loaded and then analysed as described in the manuscript and suplementary information of the associated study.
  
  * In the folder "**Auxiliary R functions**"
    - **o2cell_light_local_interpol.R**:
        This function aims at delivering a 3D matrix of light distribution in the Qubit O2 cell set up according to the optical density of the algae solution tested and sofware controlled inicident light intensity.
    - **o2cell_PO2_theo_calculation_Pm_massic_lightInhibition_modelModif.R**:
        This script aims at computing the O2 productivity in the O2 cell on single experiment knowing the kinetic parameters based on Equation 1 of the manuscript.
    - **o2cell_prod_fun_kl_corrected.R**:
        The function o2_cell_fun_Kl_corrected is defined to compute algae productivity and respiration rate as measured with the Qubit O2 cell set-up. This calculation takes into account the reaeration rate measured in the O2 cell.
    - **pc_interp2.R**:
         This function aims at accomplishing bidimensional interpolation.


The data associated needed for the code to run has been registered at the following URL: _https://doi.org/10.6084/m9.figshare.c.7746242.v1_

In order to run, the present folder can be fully downloaded. The files available at the given DOI must also be downloaded and organised as follow to run the analysis:
  - A folder '**Data files**' must be created in the root folder '**Code and data**' available in the present github repository containing the R files.
      * The file "**Results_PhoCAssPaNS.xlsx**" must be downloaded from the collection cited above and copied in the newly created folder 'Data files'. This file is available at the following [DOI:  https://doi.org/10.6084/m9.figshare.28714259]
  - A folder '**Exponential PI curve_data**' must be created in the folder '**Data files**'.
      * The files available in the item '_Results of the determination of the PI curve for Chlorella vulgaris NIES 227 in exponential growth phase_' [DOI:  https://doi.org/10.6084/m9.figshare.28714355] must be downloaded and placed in the folder 'Exponential PI curve_data'.
   - A folder '**O2Cell_data**' must finally be created in the folder '**Data files**'.
      * The files available in the item  '_Results associated to the specific oxymetric experiments performed for the publication 'Scalable modelling of photosynthetic carbon assimilation and partition during nitrogen starvation for the microalga Chlorella vulgaris NIES 227_' [DOI:  https://doi.org/10.6084/m9.figshare.28714391] must be dowloaded and placed in the folder '**O2Cell_data**'.
    
All the codes available should now be running normally.
