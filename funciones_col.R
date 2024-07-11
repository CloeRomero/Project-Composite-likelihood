# Funcion para clasificar los datos en bloques columnas
# según lo propuesto por Mukerjee (2023) en su paper
# Improving upon the effective sample size based on Godambe information for block likelihood inference
class_block_col = function(m1,m2,b1,b2,coord){
  u1= 1:m1
  u2= 1:m2
  j1= 1:b1
  j2= 1:b2
  
  # Para crear Bu1
  Bu1 =list()
  for(u in u1){
    aux=c()
    for(j in j1){
      aux= c(aux ,u+(j-1)*m1)
    }
    Bu1[[length(Bu1)+1]] = aux
  }
  
  # Para crear Bu2
  Bu2 =list()
  for(u in u2){
    aux=c()
    for(j in j2){
      aux= c(aux ,u+(j-1)*m2)
    }
    Bu2[[length(Bu2)+1]] = aux
  }
  
  # Crear el conjunto de bloques
  Bu1u2=list()
  for (i in Bu1){
    for (j in Bu2){
      Bu1u2[[length(Bu1u2)+1]] = expand.grid(i,j)
    }
  }
  
  orden = c()
  for (row_coord in 1:dim(coord)[1]){
    for(block in 1:length(Bu1u2)){
      for (row in 1:dim(Bu1u2[[block]])[1]){
        flag = all(as.numeric(Bu1u2[[block]][row,]) == as.numeric(coord[row_coord,]))
        if (flag==TRUE){
          orden = c(orden,block)
          break
        }
      }
    }
  }
  return(orden)
}


#Funcion que calcula la matriz de covarianza de Matern
# el parámetro derivative indica la derivada con respecto a qué parametro de la funcion de Matern
# 0: sin derivar, 1: derivada respecto a sigmaq, 2: derivada respecto a phi, 3: derivada respecto a tauq
calc_matrix_cov_matern = function(h, p_sigmaq, p_phi, p_tauq, derivative){
  
  if (derivative == 0){
    cov_m = p_sigmaq*exp(-h/p_phi)
  }
  if (derivative == 1 ){
    #cov_m = (1+h/p_phi)*exp(-h/p_phi)
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



# Algoritmo Fisher
Fisher_algorithm = function(theta){
  
  
  R = length(theta)
  
  print(theta)
  for (i in 1:8){
    u = array(0,R)
    H = matrix(0,R,R)
    
    for (k in 1:(M-1)){
      for (l in neighbour_list[[k]]){
        sigma_kl = calc_matrix_cov_matern(h,theta[1],theta[2],tauq,0)
        L_kl = t(chol(sigma_kl))
        Q_kl = t(solve(L_kl))%*%solve(L_kl)
        #Q_kl = try(solve(sigma_kl), stop("Sigma_kl no es invertible"))
        #tryCatch({Q_kl = solve(sigma_kl)}, error = function(ex) ({print("No se pue invertir")}))
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