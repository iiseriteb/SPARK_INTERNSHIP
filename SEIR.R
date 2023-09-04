library(deSolve)
library(shinySIR)
library(devtools)
library(shiny)

sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS <- -beta*I*S/(S+I+R) + mu*(S+I+R) + omega*R - mu*S   #susceptible
    dE <-  beta*I*S/(S+I+R) - sigma*E - mu*E   #exposed
    dI <-  sigma*E+ - gamma*I- mu*I    #infected
    dR <-  gamma*I - omega*R - mu*R     #recovery 
    return(list(c( dS, dE, dI, dR)))
  })
}


parameters_values <- c(
  mu = 1/(76*365) ,      #birth rate
  omega = 1/365,     #immunity loss
  beta  = 0.21, # infectious contact rate (/person/day)
  sigma = 1/3,    #latent period
  gamma = 1/14,          #incubation period
  delta = 0.014   #death rate
)

initial_values <- c(
  S = 999,  # number of susceptibles at time = 0
  E = 1,    #number of exposed at time =0 
  I = 0,  # number of infectious at time = 0
  R = 0   # number of recovered (and immune) at time = 0
)

time_values <- seq(0,300 , by=50 ) # days

sir_values_1 <- ode(
  y = initial_values,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)

sir_values_1

run_shiny(model = "SEIRS (w/out demography)", 
          neweqns = sir_equations ,
          ics = initial_values ,
          parm0 = parameters_values ,
          parm_names = c("birth rate", "immunity loss" , "contact rate", 
                         "latent period" ,"recovery rate" ,"death rate"),
          parm_min = c(mu=0 ,  omega = 1/365, beta = 0, sigma=1/10 , 
                       gamma=1/21, delta=0),
          parm_max = c(mu=1, omega=1 , beta = 1, sigma=1, 
                       gamma=1/7, delta=1))

with(sir_values_1, {
  # plotting the time series of susceptibles:
  plot(time, S, type = "l", col = "blue",
       xlab = "time (days)", ylab = "number of people")
  # adding the time series of infectious:
  lines(time, I, col = "red")
  # adding the time series of recovered:
  lines(time, R, col = "green")
})

# adding a legend:
legend("right", c("susceptibles", "infectious", "recovered"),
       col = c("blue", "red", "green"), lty = 1, bty = "n")

shinySIR::run_shiny("SIRS")
