


# HEX value for Kirby template code 
library(scales)
# chosing floating point scale to find the color RGB code from template 
# then using rgb function to convert to RGB code in R 

#using show the color function to see whether color is correct  
x_color <- c(rgb(0.427, 0.435, 0.447), 
             rgb(1, 0.862, 0.01),
             rgb(0.730, 0.594, 0.394), 
             rgb(0.728, 0.739, 0.746))


show_col(x_color)
