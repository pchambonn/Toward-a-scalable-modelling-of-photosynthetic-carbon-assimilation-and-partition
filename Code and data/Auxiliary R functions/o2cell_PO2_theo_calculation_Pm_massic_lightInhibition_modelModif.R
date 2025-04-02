o2cell_PO2_theo_calculation_Pm_massic_lightInhibition_modelModif <- function(I_software,Pm,K1,K2,lambda,n_mesh,OD_880,OD_683,X_algae,path){ 
  
  #Scope ------- 
  # This script aims at computing the O2 productivity in the O2 cell on single experiment knowing the kinetic parameters. 
  # The initial objective of this function is to be used for fitting of Pm and K by minimization of error functions depending on this parameters.
  # The function was modified from the model of Béchet et al 2015 to make the local light field expressed simga*X*dV*I instead of sigma*I.
  
  #Inputs----
  # I_software: light intensity as given by the software
  # Pm: Max cell specific O2 productivity (g O2/g/s)
  # K1: light saturation constant for O2 production (mumol/s/L)
  # K2 : light inhibition constant for O2 production (mumol/s/L)
  # lambda: specific respiration coefficient in the O2 cell based on same-cycle measurement (g O2/g/s) 
    # /!\ lambda is the positive coefficient for a negative contribution of respiration
  # n_mesh: number of subdivision of the O2 cell (per dimension) to use
  # OD_880: absorbance of broth at 880 nm
  # OD_683: absorbance of broth at 683 nm
  # X_algae: algae dryweight concentration (g/L)
  # path: string enabling to get the right path to the file of the function o2cell_light_local_interpol
  
  # Outputs-----
  # PO2_theo: Theoretical rate of O2 production in the O2 cell (g O2/s)
  
  # Libray and functions import----
  library(readxl)
  
  # Source functions
  
  debugSource(paste(path,'o2cell_light_local_interpol.R',sep = ""), echo=TRUE)
  
  ## Calculation----
  # Determination of light distribution and geometrical parameters:
  if (is.na(I_software) == 1){
    P_O2 <- NA
  }else{
    I_cell_table = o2cell_light_local_interpol(n = n_mesh, I_software = I_software , OD_880 = OD_880, path = path)
    
    
    
    # The coordinates (spherical) in the O2 cell are retrieved
    I_cell_model = I_cell_table$I_cell
    
    rcoordinates = I_cell_table$rcoordinates
    tetacoordinates = I_cell_table$tetacoordinates
    zcoordinates = I_cell_table$zcoordinates
    
    # We now create an array of the elementary volume for the definition used in light distribution calculation. NB: This elementary volume array as one less unit of length for the 3 coordinates.
    # An extra line of NA values for each dimension was voluntarily left in elementary volume mesh as it simplifies the restriction of I_cell on the inner mesh
    
    dV <- array(NA,dim = c(n_mesh,n_mesh,n_mesh + 1))
    
    for (i in 1:(dim(dV)[1])){
      for (j in 1:(dim(dV)[2])){
        for (k in 1:(dim(dV)[3]-1)){
          dV[i,j,k] <- rcoordinates[i]*(rcoordinates[i+1]-rcoordinates[i])*(tetacoordinates[j+1]-tetacoordinates[j])*(zcoordinates[k+1]-zcoordinates[k]) 
        }
      }
    }
    
    for (i in 1:(dim(dV)[1])){
      for (j in 1:(dim(dV)[2])){
        dV[i,j,n_mesh+1] <- rcoordinates[i]*(rcoordinates[i+1]-rcoordinates[i])*(tetacoordinates[j+1]-tetacoordinates[j])*(zcoordinates[n_mesh+1]-zcoordinates[n_mesh])
      }
    }
    
    rcoordinates <- rcoordinates[2:length(rcoordinates)]
    tetacoordinates <- tetacoordinates[1:(length(tetacoordinates)-1)]  
    
    V = sum(dV,na.rm = TRUE)
    
    
    # Determination of Pm
    sigma_683 <- 100*log(10)*OD_683/(X_algae*10^3) # m2/g
    P_O2 <- sum((Pm*(sigma_683*(X_algae)*I_cell_model)/(K1 + sigma_683*(X_algae)*I_cell_model + (1/K2)*(sigma_683*(X_algae)*I_cell_model)^2) -  lambda)*X_algae*dV*10^3,na.rm = TRUE) # The 10^3 factor converts the volume in m3 into L
  # P_O2 is determined in g O2/s
  }
  
  return(as.numeric(P_O2))
}