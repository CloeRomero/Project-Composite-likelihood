library("geoR")
library("schoolmath")

source("funciones_col.R");


# Configuracion de la grilla
k=1          # factor para modificar el tamaño de la grilla 
n1 = 16*k;   # tamaño grilla horizontal
n2 = 8*k;    # tamaño grilla vertical
n = n1*n2;   

# Configuración de los bloques
# Tiene que calzar que m1*m2 = cantidad de bloques
# b1*b2 = tamaño de los bloques
# m1*b1 = n1 y que m2*b2=n2
m1=2
b1=8*k
m2=1
b2=8*k
M=m1*m2

#Parámetros del modelo
sigmaq=1           # sill
phi=3              # rango
#par= c(sigmaq,phi);

########### CREACION PAPER PARA GUARDAR ESTIMACIONES #########################
columns = c("sigma_real", "phi_real", "sigmaq_inicial", "phi_inicial", "sigma_est", "phi_est", "tiempo")
df = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df) = columns

############# LECTURA CAMPOS GAUSSIANOS ##############################
field = read.csv("field_128_2x1_theta1-3_500r.csv")
coord = cbind(field$V1,field$V2)
######################################################################

########### CODIGO PRUEBA ##################################
# Lista de los vecinos por bloque
neighbour_list = as.list(c(2:(M)))

# Lista con la clasificación de los datos por bloque
class = class_block_col(m1,m2,b1,b2,coord)

# Coordenadas con clasificacion por bloques
coord_class = as.data.frame(cbind(coord,class))

# Matriz de distancias
h = as.matrix(dist(rbind(coord_class[coord_class$class == 1, ][,][1:2],coord_class[coord_class$class == 2, ][,][1:2]), diag=TRUE, upper=TRUE))

t_i = Sys.time()
for (i in 1:(ncol(field)-3)){
  
  data = field[,i+3]     # Lee los datos Y de cada campo
  
  ########### CODIGO PRUEBA ##################################
  
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
t_f = Sys.time() -t_i
t_f

# Crea archivo para guardar estimaciones
write.csv(df, file=paste("RF",n,"-2x1-theta1-",phi,"-iter500_C.csv",sep=""))

####################
#####################################$##########