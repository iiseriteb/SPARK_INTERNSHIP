


library(deSolve)
library(shinySIR)
library(devtools)
library(shiny)
install_github("SineadMorris/shinySIR" , force = TRUE)



sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS <- -((beta*(A+I)*S) - (alpha*E*S))/(S+A+I+R) + mu*(S+A+I+R) + omega*(R+A)    #susceptible
    dE <-  ((beta*(A+I)*S) + (alpha*E*S))/(S+A+I+R) + beta*I*R/(S+A+I+R) -(1/k)*E -lambda*E   #exposed
    dA <-  (1/k)*E - omega*A                   #asymptomatic
    dI <-  lambda*E - rho*I - theta*I  #infected but not yet tested
    dQ <-  theta*I - rho*Q - gamma*Q    #quarantine
    dH <-  rho*(I+Q) - gamma_h*H         #hospitalization
    dR <-  gamma_h*H + gamma*Q- omega*R - beta*I*R/(S+A+I+R)      #recovery 
    dD <-  delta*(I+Q+H)                            #death
    return(list(c( dS, dE, dA, dI,dQ , dH, dR , dD)))
  })
}

parameters_values <- c(
  mu = 0.018 ,     #birth rate
  p = 0.55 ,         #vaccinated
  omega = 1/240,     #immunity loss
  beta  = 0.195, # infectious contact rate (/person/day)
  alpha = 0.135,  #exposed contact rate
  lambda = 1/3,    #latent period
  gamma = 0.075,    # recovery rate (/day)
  gamma_h = 0.055,    #recovery rate from hospital
  k = 7,          #incubation period
  theta = 0.25,     #quarantine rate
  rho = 0.70,    #hospitalization rate
  delta = 0.014,   #death rate
  phi = 10       #foreigner
  )

initial_values <- c(
  #N = sum(as.numeric("S","E","A","I","R")),
  S = 100,  # number of susceptibles at time = 0
  E = 0,    #number of exposed at time =0 
  A = 0,    #number of asymptomatic at time=0
  I = 5,  # number of infectious at time = 0
  Q = 0,   #number of quarantization at time=0
  H = 0,   #number of hospitalization at time=0
  R = 0,   # number of recovered (and immune) at time = 0
  D = 0    #number of death at time=0
  )

time_values <- seq(0,300 , by=30 ) # days

sir_values_1 <- ode(
  y = initial_values,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)

shinySIR::seirs.app

run_shiny(model = "SEIRS", 
          neweqns = sir_equations ,
          ics = initial_values ,
          tstart = 0,
          timestep = 1, tmax = 300,
          parm0 = parameters_values ,
          parm_names = c("birth rate", "vaccinated" , "immunity loss" , "contact rate","contact rate from exposed", 
                         "latent period" , "recovery rate", "recovery rate from hospital" ,
                         "incubation period","testing efficiency" , "hospitalization rate" , "death rate" , "foreigners"),
          parm_min = c(mu=0 , p=0 , omega = 1/365, beta = 0, alpha = 0, lambda=1/10 , 
                       gamma = 0 , gamma_h=0, k=7, theta=0 , rho=0 , delta=0, phi=0),
          parm_max = c(mu=round(1/365 ,3), p=1, omega=1 , beta = 1, alpha=1, lambda=1, 
                       gamma = 1 , gamma_h=1 ,k=21, theta=1 , rho=1 , delta=1, phi=50),
          sigfigs = 4, showtable = TRUE, linesize = 1.2, textsize = 14 ,
          legend_title = "Compartment",
slider_steps = NULL, values = NULL,)
          

R0 <- as.numeric( beta^2 *lambda / (omega*k *(rho + theta) *(lambda + 1/k)^2))


