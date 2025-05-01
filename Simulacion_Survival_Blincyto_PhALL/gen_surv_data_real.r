
# Funcion para generar tiempos de supervivencia y censuras con 
# la aleatorizacion hecha anteriormente usando la distribucion Weibull
# Este codigo ya tiene en cuenta la posibilidad de que hayan mas de dos brazos (o uno solo)

gen_surv_data_real <- function(pat, parameters, 
                          shape_parameter,
                          censor, Tmax){
  
  # Identificamos los diferentes brazos que hay
  
  arms <- sort(unique(pat$arm))
  
  # Inicializamos vectores
  
  event_times <- numeric(length(pat$id))
  censor_times <- numeric(length(pat$id))
  
  # Bucle por nº tratamientos
  
  for (arm in arms) {
    
    # Cogemos el nº de pacientes del brazo en el que estamos
    
    n <- sum(pat$arm == arm)
    
    # Generar eventos con esta distribucion
    
    event_times[pat$arm == arm] <- rweibull(n, shape_parameter[arm], parameters[arm])
    
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
  
  
  return(pat)
  
}