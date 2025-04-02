o2cell_light_local_interpol <- function(n,I_software, OD_880, path){
  
  source(paste(path,'/pc_interp2.R',sep = ""), echo=TRUE)

  # This function aims at delivering a 3D matrix of light distribution in the O2 cell set up according
  # to the optical density of the algae solution tested and sofware controlled inicident light intensity.
  # The calculation is developed from measurements performed on 28/09/21 which enabled to "model" light on
  # the sides and on the central axis of the O2 cell. Light intensity in the rest of the volume is computed
  # from linear interpolation between the horizontal projections on the central axis and the side of the O2 cell.
  # The calculation is made in cylindrical coordinates.
  
  ## INPUTS
  # n: number of subdivision of the volume in each dimension
  # I_sofware: light intensity as delivered according to the O2 cell software (mumol/m2/s)
  # OD_880: optical density of the algae solution measured at 880 nm in 1cm large cells.
    ## NB: this measurement was chosen as it was found to enable to predict the extinction of the light 
    ## measured in the O2 cell with a fair accuracy. However, algae of varying health may induce an error
    ## in the total light absorbed based on the relations developed (see 21_09_27_light_distribution.xls )
  # path : path to the file containing pc_interp2 function.
  
    
  ## OUTPUTS
    
    
  ## Parameters definition--------------------
  R = 1.27/2*10^(-2)
  S_top = pi*R^2
  h = 5.51*10^(-6)/S_top
  kd = 1.5 # empirical "light extinction coefficient" for inside and exit measurement points  
    
  ## Light intensity pre-definition--------------------
  
  # I0_meas_100 is the array of the light intensity measured 
    # Matrices of sample points
    # Dimension 1: sigma 0:  12,3,12'; 11,2,11'; 10,1,10'
    # Dimension 2: sigma_in: 12,6,12'; 11,5,11'; 10,4,10'
    # Dimension 3: sigma_exit: 12,9,12'; 11,8,11'; 10,7,10'.
  
  I0_meas_100 <- array(data =c(39,51,39, 91,115,91, 6,6,6, 
                               39,36,39, 91,105,91, 6,7,6, 
                               39,29,39, 91,91,91, 6,9,6),
                       dim = c(3,3,3))
   
  # The matrix of known intensity is now corrected for the incident light intensity delivered based on 
  # software indication and the light extinction coefficient at 880 nm of the algae solution.
  
  # NB:
  I0_corrected <- I_software/101.1*
                  array(data = c(max(c(0,36.0-10.8*OD_880)) , max(c(0,55.1-17.7*OD_880)) , max(c(0,36.0-10.8*OD_880)), # Keep adding the max to all the linear projections
                                 max(c(0,90.6-44.0*OD_880)) , max(c(0,118-13.87*OD_880)) , max(c(0,90.6-44.0*OD_880)),
                                 max(c(0,5.72-3.46*OD_880)) , max(c(0,5.32-1.66*OD_880)) , max(c(0,5.72-3.46*OD_880)),
                                 
                                 max(c(0,36.0-10.8*OD_880)) , I0_meas_100[2,1,2]*exp(-kd*OD_880) , max(c(0,36.0-10.8*OD_880)),
                                 max(c(0,90.6-44.0*OD_880)) , I0_meas_100[2,2,2]*exp(-kd*OD_880) , max(c(0,90.6-44.0*OD_880)),
                                 max(c(0,5.72-3.46*OD_880)) , I0_meas_100[2,3,2]*exp(-kd*OD_880) , max(c(0,5.72-3.46*OD_880)),
                                 
                                 max(c(0,36.0-10.8*OD_880)) , I0_meas_100[2,1,3]*exp(-kd*OD_880) , max(c(0,36.0-10.8*OD_880)),
                                 max(c(0,90.6-44.0*OD_880)) , I0_meas_100[2,2,3]*exp(-kd*OD_880) , max(c(0,90.6-44.0*OD_880)),
                                 max(c(0,5.72-3.46*OD_880)) , I0_meas_100[2,3,3]*exp(-kd*OD_880) , max(c(0,5.72-3.46*OD_880))),
                        dim = c(3,3,3))
  
  
  ## Interpolation of light intensity in the three plans where value is known-----
  
  # SO
  
  teta_sample_side_S0 = c(pi , 3*pi/2 , 2*pi) # angluar coordinates of the frontal face
  z_sample_side_S0= c(0 , h/3 , h) # vertical coordinates of the frontal face
  
  nteta <- n
  nz <- n
  
  interpolated_I_S0 <- pc_interp2(teta_sample_side_S0 , z_sample_side_S0 ,  I0_corrected[,,1] , nteta , nz)
    teta_mod_side_S0 <- rev(interpolated_I_S0$x) # the rev() is introduced because on this face the data in the matrix is in the reverse order of angles (written from 2*pi to pi)
    z_mod_side_S0 <- interpolated_I_S0$y
    I0_mod_S0 <- interpolated_I_S0$z
    
  # S1
    
  teta_sample_side_S1 = c(0 , pi/2 , pi) # angluar coordinates of the exit face
  z_sample_side_S1= c(0 , h/3 , h) # vertical coordinates of the exit face
  
  interpolated_I_S1 <- pc_interp2(teta_sample_side_S1 , z_sample_side_S1 ,  I0_corrected[,,3] , nteta , nz)
  teta_mod_side_S1 <- interpolated_I_S1$x
  z_mod_side_S1 <- interpolated_I_S1$y
  I0_mod_S1 <- interpolated_I_S1$z
  
  
  # Central axis
  z_sample_side_axis= c(0 , h/3 , h) # vertical coordinates of the central axis
  interpolated_I_axis = approx(z_sample_side_axis,I0_corrected[2,,2], method = "linear",n = nz)
  z_mod_axis <- interpolated_I_axis$x
  I0_mod_axis <- interpolated_I_axis$y
  
  
  ## Mapping of O2 cell inner coordinates------------------------
  
  nr = n
  nteta = n
  nz = n
  
  rcoordinates = seq(0,R,R/nr)
  rcoordinates <- rcoordinates[2:length(rcoordinates)]
  tetacoordinates = seq(0,2*pi,2*pi/nteta)
  tetacoordinates <- tetacoordinates[1:(length(tetacoordinates)-1)]
  zcoordinates = seq(0,h,h/nz)
  
  nr = length(rcoordinates)
  nteta = length(tetacoordinates)
  nz = length(zcoordinates)
  
  
  ## Determination of O2 cell light distribution----------------
  
  I_cell <-  array(,dim = c(nr,nteta,nz))
  
  for (ir in 1:nr){
    for (iteta in 1:nteta){
      for (iz in 1:nz){
        I_axis <- I0_mod_axis[which.min(abs(z_mod_axis - zcoordinates[iz]))]
        if (tetacoordinates[iteta] <= pi){
          I_side <- I0_mod_S1[which.min(abs(teta_mod_side_S1 - tetacoordinates[iteta])),
                              which.min(abs(z_mod_side_S1 - zcoordinates[iz]))]
        }else{
          I_side <- I0_mod_S0[which.min(abs(teta_mod_side_S0 - tetacoordinates[iteta])),
                              which.min(abs(z_mod_side_S0 - zcoordinates[iz]))]
        }
        I_cell[ir,iteta,iz] = approx(c(0,R),c(I_axis,I_side),xout = rcoordinates[ir])$y
        
        
      }
    }
  }
  
  outputs <- list(I_cell = I_cell,
                  rcoordinates = rcoordinates , tetacoordinates = tetacoordinates , zcoordinates = zcoordinates,
                  I_0_mapped = I0_corrected)

}