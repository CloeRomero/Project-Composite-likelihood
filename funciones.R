#install.packages("geoR")
library(geoR)
#install.packages("Matrix")
library(Matrix)
#install.packages("schoolmath")
library("schoolmath")

####################################################################
###################### FUNCIONES ###################################
####################################################################


# Funcion que clasifica los datos por bloque con el enfoque fila
class_blocks = function(coord, nb_x, nb_y, x_size, y_size){
  class = c()
  block_x_size = x_size/nb_x
  block_y_size = y_size/nb_y
  
  for(i in seq(n)){
    x_pos = (as.numeric(coord[i,][1])-1)%/%block_x_size
    y_pos = (as.numeric(coord[i,][2])-1)%/%block_y_size
    
    if (x_pos >= nb_x){
      x_pos = x_pos-1
    }
    if (y_pos >= nb_y){
      y_pos = y_pos-1
    }
    
    block_number = 1 + x_pos + nb_x*y_pos
    #print(block_number)
    class = c(class,block_number)
  }
  return(class)
}

# Funcion que entrega los vecinos hacia delante de cada bloque
neighbour_map = function(nb_x,nb_y){
  neighbourmap = list()
  for (i in seq(nb_x*nb_y-1)){
    neighbours = c()
    if (i < nb_x*(nb_y-1)+1){
      neighbours = c(neighbours,i + nb_x)
      
      if ((i-1)%%nb_x != 0){ # si no está a la izquierda
        neighbours = c(neighbours,i + nb_x -1)
      }
      if (i%%nb_x != 0){ # si no está a la derecha
        neighbours = c(neighbours,i + nb_x + 1, i + 1)
      }
    } else{
      if (i%%nb_x != 0){ # si no está a la derecha
        neighbours = c(neighbours,i + 1)
      }
    }
    
    neighbourmap = c(neighbourmap, list(neighbours))
  }
  
  return(neighbourmap)
}

#Funcion que calcula la matriz de covarianza de Matern
# el parámetro derivative indica la derivada con respecto a qué parametro de la funcion de Matern
# 0: sin derivar, 1: derivada respecto a sigmaq, 2: derivada respecto a phi, 3: derivada respecto a tauq
calc_matrix_cov_matern = function(h, p_sigmaq, p_phi, p_tauq, derivative){
  
  if (derivative == 0){
    cov_m = p_sigmaq*exp(-h/p_phi)
  }
  if (derivative == 1 ){
    cov_m = exp(-h/p_phi)
  }
  if (derivative == 2){
    cov_m = p_sigmaq*h*exp(-h/p_phi)/(p_phi^2)
  }
  if (derivative == 3){
    cov_m = ifelse(h==0,1,0)
  }
  
  return(cov_m)
}


# Funcion que calcula las matrices de distancia de cada par de bloques vecinos hacia adelante
# Esta funcion solo es valida para configuraciones de bloque donde todos los bloques sean de igual tamaño y forma
# En la posición 1 de la lista queda la matriz de distancia entre bloques adyacentes uno a la derecha del otro
# En la posición 2*nb_x de la lista queda la matriz de distancia entre un bloque y su bloque adayacente que está arriba
# En la posición 2*nb_x-1 de la lista queda la matriz de distancia entre un bloque y su bloque adyacente diagonal arriba a la izquierda
# En la posición 2*nb_x+1 de la lista queda la matriz de distancia entre un bloque y su bloque adyacente diagonal arriba a la derecha
matrices_distancia = function(nb_x, nb_y){
  matrices_h = as.list(rep(0, 2*nb_x+1))
  if (nb_y == 1){
    m_dist = as.matrix(dist(rbind(coord_class[coord_class$class == 1, ][,][1:2],coord_class[coord_class$class == 2, ][,][1:2]), diag=TRUE, upper=TRUE));
    matrices_h[[ 1 ]] = m_dist
  }else{
    if(nb_x == 2){
      #bloque a la derecha
      matrices_h[[1]] = as.matrix(dist(rbind(coord_class[coord_class$class == 1, ][,][1:2],coord_class[coord_class$class == 2, ][,][1:2]), diag=TRUE, upper=TRUE));
      #bloque arriba izquierda
      matrices_h[[2*nb_x - 1]] = as.matrix(dist(rbind(coord_class[coord_class$class == 2, ][,][1:2],coord_class[coord_class$class == 3, ][,][1:2]), diag=TRUE, upper=TRUE));
      #bloque arriba
      matrices_h[[2*nb_x]] = as.matrix(dist(rbind(coord_class[coord_class$class == 1, ][,][1:2],coord_class[coord_class$class == 3, ][,][1:2]), diag=TRUE, upper=TRUE));
      #bloque arriba derecha
      matrices_h[[2*nb_x + 1]] = as.matrix(dist(rbind(coord_class[coord_class$class == 1, ][,][1:2],coord_class[coord_class$class == 4, ][,][1:2]), diag=TRUE, upper=TRUE));
    }else{
      #bloque a la derecha
      matrices_h[[1]] = as.matrix(dist(rbind(coord_class[coord_class$class == 1, ][,][1:2],coord_class[coord_class$class == 2, ][,][1:2]), diag=TRUE, upper=TRUE));
      #bloque arriba izquierda
      matrices_h[[2*nb_x - 1]] = as.matrix(dist(rbind(coord_class[coord_class$class == 2, ][,][1:2],coord_class[coord_class$class == (2+nb_x-1), ][,][1:2]), diag=TRUE, upper=TRUE));
      #bloque arriba
      matrices_h[[2*nb_x]] = as.matrix(dist(rbind(coord_class[coord_class$class == 2, ][,][1:2],coord_class[coord_class$class == (2+nb_x), ][,][1:2]), diag=TRUE, upper=TRUE));
      #bloque arriba derecha
      matrices_h[[2*nb_x + 1]] = as.matrix(dist(rbind(coord_class[coord_class$class == 2, ][,][1:2],coord_class[coord_class$class == (2+nb_x+1), ][,][1:2]), diag=TRUE, upper=TRUE));
    }
  }
  return(matrices_h)
}


# Funcion que calcula las inversas de las matrices de covarianza a partir de las matrices de distancia
# calculadas en la funcion "matrices_distancia". La posición de las matrices inversas en la lista resultante
# sigue la misma lógica que la de las matrices de distancia
calc_matrices_Q_kl = function(nb_x, nb_y, theta, m_dist){
  matrices_Q_kl = as.list(rep(0, 2*nb_x + 1))
  
  if (nb_y == 1){
    sigma_kl = calc_matrix_cov_matern(m_dist[[1]],theta[1],theta[2],tauq,0)
    L_kl = t(chol(sigma_kl))
    matrices_Q_kl[[1]] = t(solve(L_kl))%*%solve(L_kl)
  }else{
    #bloque a la derecha
    sigma_kl = calc_matrix_cov_matern(m_dist[[1]],theta[1],theta[2],tauq,0)
    L_kl = t(chol(sigma_kl))
    matrices_Q_kl[[1]] = t(solve(L_kl))%*%solve(L_kl)
    
    #bloque arriba izquierda
    sigma_kl = calc_matrix_cov_matern(m_dist[[2*nb_x - 1]],theta[1],theta[2],tauq,0)
    L_kl = t(chol(sigma_kl))
    matrices_Q_kl[[2*nb_x - 1]] = t(solve(L_kl))%*%solve(L_kl)
    
    #bloque arriba
    sigma_kl = calc_matrix_cov_matern(m_dist[[2*nb_x]],theta[1],theta[2],tauq,0)
    L_kl = t(chol(sigma_kl))
    matrices_Q_kl[[2*nb_x]] = t(solve(L_kl))%*%solve(L_kl)
    
    #bloque arriba derecha
    sigma_kl = calc_matrix_cov_matern(m_dist[[2*nb_x + 1]],theta[1],theta[2],tauq,0)
    L_kl = t(chol(sigma_kl))
    matrices_Q_kl[[2*nb_x + 1]] = t(solve(L_kl))%*%solve(L_kl)
    
  }
  return(matrices_Q_kl)
  
}


# Algoritmo Fisher

Fisher_algorithm = function(theta){
  
  R = length(theta)
  
  print(theta)
  for (i in 1:8){
    u = array(0,R)
    H = matrix(0,R,R)
    
    matrices_Q_kl = calc_matrices_Q_kl(nb_x, nb_y, theta, m_dist)
    
    
    for (k in 1:(M-1)){
      for (l in neighbour_list[[k]]){
        
        h = m_dist[[(l-k)+((l-1)%/%nb_x -(k-1)%/%nb_x)*nb_x]]
        
        Q_kl = matrices_Q_kl[[(l-k)+((l-1)%/%nb_x -(k-1)%/%nb_x)*nb_x]]
        #Q_kl = try(solve(sigma_kl), stop("Sigma_kl no es invertible"))
        #tryCatch({Q_kl = solve(sigma_kl)}, error = function(ex) ({print("No se puede invertir")}))
        q_kl = Q_kl%*%data.matrix(rbind(datos_class[datos_class$class == k,][1],datos_class[coord_class$class == l,][1]))
        
        for (r in 1:R){
          dsigma_dr = calc_matrix_cov_matern(h,theta[1],theta[2],tauq,r)
          W_klr = Q_kl%*%dsigma_dr
          #u= u + (1/2)*as.numeric(-sum(diag(W_klr)) + t(q_kl)%*%dsigma_dr%*%q_kl) # Exponencial
          u[r] = u[r] + (1/2)*as.numeric(-sum(diag(W_klr)) + t(q_kl)%*%dsigma_dr%*%q_kl) # Matern
          
          for(s in r:R){
            dsigma_ds = calc_matrix_cov_matern(h,theta[1],theta[2],tauq,s)
            W_kls = Q_kl%*%dsigma_ds
            #H=H+(1/2)*sum(diag(W_klr%*%W_kls)) #Exponencial
            H[r,s] = H[s,r] = H[r,s] + (1/2)*sum(diag(W_klr%*%W_kls)) # Matern
          }
          
        }
      }
    }
    
    theta = theta + solve(H)%*%u
    #print(theta)
  }
  return(theta)
}