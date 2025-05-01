
# Funcion para generar tiempos de supervivencia y censuras con 
# la aleatorizacion hecha anteriormente usando la distribucion Weibull
# Este codigo ya tiene en cuenta la posibilidad de que hayan mas de dos brazos (o uno solo)

gen_surv_data <- function(pat, parameters, 
                          shape_parameter,
                          censor, Tmax){
  
  # Identificamos los diferentes brazos que hay
  
  arms <- sort(unique(pat$arm))
  
  # Inicializamos vectores
  
  event_times <- numeric(length(pat$id))
  censor_times <- numeric(length(pat$id))
  
  ## SECCIÓN SÓLO PARA SIMULAR LOS DATOS REALES CON DOS DIFERENTES DIST. WEIBULL ##
  
  # COMENTAR LOS SIGUIENTE PARA SIMULACIONES SAP YA QUE IGNORA LOS INPUT DE SHAPE Y EFECTOS DE LA FUNCIÓN PRINCIPAL #
  
  # Para sacar estos valores hay que hacer lo del principio del informe con el IPD y ajustando el modelo:
  
  # Valor del shape para cada uno de los brazos
  
  # shape_parameter <- c(1.369615, 0.8690831)
  # 
  # # Valor del scale para cada uno de los brazos
  # 
  # parameters <- c(6.761968, 14.88501)
  # 
  #################################################################################
  
  # Bucle por nº tratamientos
  
  for (arm in arms) {
    
    # Cogemos el nº de pacientes del brazo en el que estamos
    
    n <- sum(pat$arm == arm)
    
    # Generar eventos con esta distribucion
    
    # Poner uno u otro en función de si simulamos datos del SAP o reales para ver si hay diferentes shapes o uno en común
    event_times[pat$arm == arm] <- rweibull(n, shape_parameter, parameters[arm]) #Datos SAP
    #event_times[pat$arm == arm] <- rweibull(n, shape_parameter[arm], parameters[arm]) # Datos reales
    
    
    # Generamos ahora tiempos de censura
    
    # Estimamos la proporcion de pacientes (1-censor[arm]) que van a experimentar un EVENTO
    # Para ello calculamos los cuantiles de los tiempos estimados de eventos pero esto con la
    # distribucion exponencial
    
    censor_times <- rexp(n, rate = 1 / (quantile(event_times[pat$arm == arm], 1 - censor[arm])))
  
    }
  
  # Aplicamos censura. En status:
  
  # TRUE = Event
  # FALSE = Censura
  
  status <- (event_times <= censor_times) & (event_times <= Tmax) 
  
  # Tiempo de evento o censura
  
  Event_Censored <- pmin(pmin(event_times, censor_times), Tmax)
  
  # Aplicamos censura administrativa (i.e., tiempos de evento que se pasan del tiempo max de estudio)
 
   status[Event_Censored == Tmax] <- FALSE
  
  ## PARA HACERLO MAS REAL HACE FALTA HACER LO DEL TIEMPO DE LLEGADAS POR PACIENTE!!
  

  # Incorporamos estos datos al df pat.
   
  pat$status <- status
  pat$time <- Event_Censored
  
  # Mirar esto para mañana yu ver que falla
  
  return(pat)
  
}