# THis function aims at accomplishing bidimensional interpolation. Knowing the values 
# of a variable z on (x,y) coordinates, this function returns the hypothetical values 
# z on coordinates from linear interpolation of x and y data.  


pc_interp2 <- function(x , y , z , nx , ny  ){
  
  # INPUTS
  
  # x: intial x coordinates (ordered values, dimension lx)
  # y: inital y coordinates (ordered values, dimension ly)
  # z: know values of the surface to interpolate (dimension lx by ly) or ly by lx???? unclear
  # nx: number of point by which to divide in between each known points of x
  # ny: number of point by which to divide in between each known points of y
  
  # OUTPUTS
  
  # x_interp: vector of interpolated data of x: nx - 1 values are created in between point x. The lenght of x_interp is (nx-1)*length(x) +1
  # y_interp: vector of interpolated data of y: ny - 1 values are created in between each points of y. The lenght of y_interp is (ny-1)*length(y) +1

  x_interp <- c(x[1])
  for (i in 1:(length(x) - 1)){
    x_interp <- c(x_interp ,  approx(x[i:(i+1)] , method = "linear" , n = nx + 1 )$y[-1])  
  }
  y_interp <- c(y[1])
  for (i in 1:(length(y) - 1)){
    y_interp <- c(y_interp ,  approx(y[i:(i+1)] , method = "linear" , n = ny + 1)$y[-1])
  }
  
  z_interp <- matrix( , nrow = length(x_interp) , ncol = length(y_interp))
  for (ix in 1:(length(x_interp)-1)){
    for (iy in 1:(length(y_interp)-1)){
      # Find framing values of x and y and respective index
      x1 <- max(x[x <= x_interp[ix]])
      index_x1 <- which(x == max(x[x <= x_interp[ix]]))
      x2 <- x[index_x1 + 1]
      index_x2 <- index_x1 + 1
      
      y1 <- max(y[y <= y_interp[iy]])
      index_y1 <- which(y == max(y[y <= y_interp[iy]]))
      y2 <- y[index_y1 + 1]
      index_y2 <- index_y1 + 1
                        
      z_interp[ix,iy] <- z[index_x1 , index_y1] + 
        (z[index_x2 , index_y1] - z[index_x1 , index_y1])*(x_interp[ix] - x1)/(x2 - x1) +
        (z[index_x1 , index_y2] - z[index_x1 , index_y1])*(y_interp[iy] - y1)/(y2 - y1) +
        (z[index_x2 , index_y2] + z[index_x1 , index_y1] - z[index_x2 , index_y1] - z[index_x1 , index_y2])*
          (x_interp[ix] - x1)/(x2 - x1)*(y_interp[iy] - y1)/(y2 - y1)
    }
  }
  # adding last row and last column
  last_col <-  c(z[1,ncol(z)])
  for (i in 1:(length(x) - 1)){
    last_col <- c(last_col ,  approx(z[i:(i+1),ncol(z)] , method = "linear" , n = nx + 1)$y[-1])
  }
  z_interp[,ncol(z_interp)] <- last_col
  
  last_row <- c(z[nrow(z),1])
  for (i in 1:(length(y) - 1)){
    last_row <- c(last_row ,  approx(z[nrow(z),i:(i+1)] , method = "linear" , n = ny + 1)$y[-1])
  }
  z_interp[nrow(z_interp),] <- last_row
    
  interpolated_data <<- list(x = x_interp, y = y_interp, z = z_interp)
}