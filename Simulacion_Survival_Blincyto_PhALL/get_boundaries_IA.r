
# Funcion para estraer los valores de p limites para declarar si es significativo el resultado o no
# Solo extraen los limites para un IA dado, de momento solo trabajamos con O'Brien-Fleming y Pocock


get_boundaries_IA <- function(IA, alpha, method_IA = "OF", test.type = test.type,
                              sample_size, n_exp_events){
  
  if (method_IA == "OF") {
    
    
    x <- gsDesign(k=length(IA)+1, sfu="OF", timing = IA,
                                 test.type = 2, alpha = alpha, 
                                 nFixSurv = sample_size, n.fix = n_exp_events)
    

    boundary <- pnorm(q=x$upper$bound, lower.tail=FALSE) 
    
    return(boundary)
    
  } 
    
  if (method_IA == "P") {
    
    x <- gsDesign(k=length(IA)+1, sfu="Pocock", timing = IA,
                  test.type = 2, alpha = alpha,
                  nFixSurv = sample_size, n.fix = n_exp_events)
    
    
    boundary <- pnorm(q=x$upper$bound, lower.tail=FALSE)
    
    return(boundary)
    
  }
  
  if (method_IA == "Bayes") {
    
    boundary <- rep(alpha, length(IA)+1)
    
    return(boundary)
    
  } 
  
}

