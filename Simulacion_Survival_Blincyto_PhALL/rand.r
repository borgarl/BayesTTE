# Funcion para aleatorizar pacientes:
# 1) CR -> Complete randomisation.
# 2) BR -> Block randomisation.


# PARAMETROS:
# sample_size -> Numero total de pacientes.
# ratio -> ratio de aleatorizacion (ej: 0.5,0.5 es 1:1)
# block_size -> Si BR=TRUE, block_size es el tamaño de cada uno de los bloques.
# CR -> Indicador binario que indica el metodo de aleatorizacion escogido.

rand <- function(sample_size, ratio, rand_type, block_size = FALSE) {
  
  # Calculo de nº de brazos de tratamiento
  n_treatments <- length(ratio)
  
  # Calculo de nº de pacientes por cada brazo de tratamiento
  # round en caso de que la division no sea entera
  n_per_treatment <- round(sample_size*ratio)
  
  
  # Creamos un vector para la asignacion de pacientes por cada brazo
  trt  <- rep(1:n_treatments, times = n_per_treatment)
  
  # Inicializar vector
  
  assign <- data.frame(matrix(trt, nrow = length(trt), ncol = 2)) 
  columns = c("id","arm") 
  colnames(assign) = columns

  
  # Antes de asignar brazos,se genera un ID para cada paciente simulado 
  # (esto quizás sea innecesario para muchas simulaciones)
  prefix <- "PAT"
  rand_num <- sample(1:length(trt), length(trt), replace = FALSE)
  id <- paste0(prefix, rand_num)
  assign$id <- id
  
  
  # En caso de que BR==TRUE se hace la aleatorizacion por bloques.
  # Esto da error si el tamaño de cada bloque no deja resto=0 en el tamaño muestral
  # Modificaremos si hace falta pero no ahora
  if (rand_type == "BR") {
    
    # Se calcula el nº de bloques en funcion del tamaño de cada bloque y el total de la muestra
    num_blocks <- ceiling(length(trt) / block_size)
    
    # Se ordenan
    blocks <- rep(1:num_blocks, each = block_size)
    
    # En caso de que salga con decimal lo anterior
    if (length(trt) %% block_size != 0) blocks <- c(blocks, rep(tail(blocks, 1), length(trt) %% block_size))
    
    # Aleatorizar pacientes en cada uno de los bloques
    rand_pat <- rep(NA, length(trt))
    
    for (i in unique(blocks)) {
      
       block_i <- blocks == i


      n_patients_per_arm <- round(sum(block_i) * ratio / sum(ratio))
      
      n_patients_per_arm[n_treatments] <- sum(block_i) - sum(n_patients_per_arm[-n_treatments])
      
      treatments <- unlist(lapply(1:n_treatments, function(arm) rep(arm, n_patients_per_arm[arm])))
      
      rand_pat[block_i] <- sample(treatments)
      
    }
    
  
    assign$arm <- rand_pat
    
  }
  
  # En caso de que sea FALSO (está por defecto en caso de no especificarlo)
  
  else {

    # Simplemente "barajamos" de manera aleatoria el vector (esto vale para ECs con cualquier n de brazos)
    assign$arm <- sample(assign$arm)
    
  }
  
  
  
  # Devuelve un vector con los pacientes aleatorizados por cada brazo 
  
  return(assign)
  
}
