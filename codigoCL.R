#install.packages("geoR")
library(geoR)
#install.packages("Matrix")
library(Matrix)
#install.packages("schoolmath")
library("schoolmath")

source("funciones.R")

# Configuracion de bloques
nb_x = 2         # numero de bloques horizontales
nb_y = 1         # numero de bloques verticales
M = nb_x*nb_y    # numero de bloques

# Configuracion campo Gaussiano
grid = "reg"       #Tipo de grilla
k=1                #factor para ajustar el tamaño 
nx = 16*k          # cantidad de observaciones horizontales
ny = 8*k           # cantidad de observaciones verticales
n = nx*ny
x_grid_size = nx   # ancho de la grilla
y_grid_size = ny   # alto de la grilla

sigmaq=1   # sill
phi=3      # rango
tauq=0

# Data frame para almacenar las estimaciones
columns = c("sigma_real", "phi_real", "sigmaq_inicial", "phi_inicial", "sigma_est", "phi_est", "tiempo")
df = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df) = columns

############# LECTURA CAMPOS GAUSSIANOS ##############################
field = read.csv("field_128_2x1_theta1-3_500r.csv")
coords = cbind(field$V1,field$V2)
######################################################################


########### CODIGO CL ##################################
# Lista de los vecinos por bloque
neighbour_list = neighbour_map(nb_x,nb_y)

# Lista con la clasificación de los datos por bloque
class = class_blocks(coords, nb_x, nb_y, nx, ny)

# Coordenadas con clasificacion por bloques
coord_class = as.data.frame(cbind(coords,class))

# Matrices de distancia
m_dist = matrices_distancia(nb_x, nb_y)

t_i = Sys.time()
for (i in 1:(ncol(field)-3)){ # cuantos campos hay almacenados en field
  
  data = field[,i+3]   # Lee los datos Y de cada campo
  
  ########### CODIGO CL ##################################
  
  # Datos Y con clasificación por bloques
  datos_class = as.data.frame(cbind(data,class))
  
  # Theta inicial
  theta_0 = c(sigmaq + round(runif(1,1,4)), phi + round(runif(1,1,4)))
  
  theta=c(0,0)
  flag = "No error"
  
  t <- proc.time()[[3]] # Inicia el cronómetro
  tryCatch({theta = Fisher_algorithm(theta_0)}, error = function(ex) ({flag="Error";print(ex)}))
  time=proc.time()[[3]]-t    # Detiene el cronómetro
  
  # Guarda las estimaciones en el dataframe
  df[nrow(df) + 1,] = c(sigmaq,phi,theta_0[1], theta_0[2],theta[1],theta[2], time)
  
  print(i)
}  

t_f = Sys.time() - t_i
t_f

# Crea archivo para guardar estimaciones
write.csv(df, file=paste("RF",n,"-",nb_x,"x",nb_y,"-theta1-",phi,"-iter500.csv",sep=""))

#################
#######################