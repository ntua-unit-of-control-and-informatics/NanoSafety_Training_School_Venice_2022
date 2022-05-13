###############################################################################################################
###  First step: install the jaqpotr package that allows connection with the jaqpot cloud platform        #####
###  and the deSolve R package for trying out the ODE model                                               #####
###############################################################################################################

devtools::install_github("euclia/jaqpotr", force = TRUE) 
install.packages("deSolve")
install.packages("ggplot2")

#########################################################################################
###  D.magna model R code. The model was published by Fan et al. (2016)             #####
###  and can be found at https://doi.org/10.1016/j.scitotenv.2016.06.197            #####
#########################################################################################


#=================
#1. User input 
#=================
material <- "S2"  # Choose one from T1,T2,H1,S1,H2,S2
C_water <- c(1,0,10, 15) # NP concentration in water, in mg/L. Experimentally
# tested values: 0.1, 1.0 & 10.0 mg/L
sampling_time <- c(0,5,10, 20) #sampling time, in hours
user_input <- list("material"=material, "C_water"=C_water, "sampling_time" = sampling_time )

#===========================================
#2. Function that creates parametric vector    
#===========================================
create.params <- function(user_input){
  with( as.list(user_input),{
    # Uptake rate constants (L/g/h)
    ku_data <- matrix(c(11.842, 14.551, 7.264,
                        10.042, 17.821, 11.491,
                        17.640, 19.427, 8.935,
                        28.153, 30.855, 14.369,
                        18.354, 15.312, 8.774, 
                        15.340, 32.942, 17.455), nrow = 6, byrow = TRUE)
    rownames(ku_data) <- c("T1","T2","H1","S1","H2","S2")
    colnames(ku_data) <- c("0.1","1","10")                      
    
    # Efflux rate constants (1/h)                     
    ke_data <- matrix(c(0.0257, 0.0136, 0.0228,
                        0.0424, 0.0443, 0.0184,
                        0.0341, 0.0342, 0.0157,
                        0.0327, 0.0289, 0.0229,
                        0.0234, 0.0174, 0.0146,
                        0.0207, 0.0419, 0.0119), nrow = 6, byrow = TRUE)                       
    rownames(ke_data) <- c("T1","T2","H1","S1","H2","S2")
    colnames(ke_data) <- c("0.1","1","10")  
    
      # Parameter averaging over the three concentrations of the experiment to
      # integrate out the impact of the concentration magnitude
      if (material=="T1"){
        ku <- mean(ku_data[1,1], ku_data[1,2], ku_data[1,3]) 
        ke <- mean(ke_data[1,1], ke_data[1,2], ke_data[1,3])
      } else if(material=="T2"){
        ku <- mean(ku_data[2,1], ku_data[2,2], ku_data[2,3])
        ke <- mean(ke_data[2,1], ke_data[2,2], ke_data[2,3])
      } else if(material=="H1"){
        ku <- mean(ku_data[3,1], ku_data[3,2], ku_data[3,3])
        ke <- mean(ke_data[3,1], ke_data[3,2], ke_data[3,3])
      } else if(material=="S1"){
        ku <- mean(ku_data[4,1], ku_data[4,2], ku_data[4,3])
        ke <- mean(ke_data[4,1], ke_data[4,2], ke_data[4,3])
      } else if(material=="H2"){
        ku <- mean(ku_data[5,1], ku_data[5,2], ku_data[5,3])
        ke <- mean(ke_data[5,1], ke_data[5,2], ke_data[5,3])
      } else if(material=="S2"){
        ku <- mean(ku_data[6,1], ku_data[6,2], ku_data[6,3])
        ke <- mean(ke_data[6,1], ke_data[6,2], ke_data[6,3])
      } else {stop("Select a valid material")}
      

    return(list("ku"=ku, "ke"=ke, "Cw"=C_water, "sampling_time" = sampling_time))
  }) 
}

#===============================================
#3. Function that creates initial values for ODEs 
#===============================================
create.inits <- function(parameters){
  with( as.list(parameters),{
    print(sampling_time)
    C <- 0; Cw <-0; 
    
    return(c("C"=C, "Cw"=Cw))
  })
}

#==================================
#4. Function that creates events      
#==================================

create.events <- function(parameters){
  with( as.list(parameters), {
    # Append a zero at the end of all concentrations to stop exposure
    Cw <- c(Cw,0)
    L0 <- length(sampling_time)
    # Add another hour to the sampling time vector
    sampling_time <- c(sampling_time, sampling_time[L0]+1)
    lCw <- length(Cw)
    ltimes <- length(sampling_time)
    if (lCw != ltimes){
      stop("The length of C_water and sampling_time must be equal.")
    }else{
      events <- list(data = data.frame(var = c(rep("Cw", lCw)),
                                       time = sampling_time, 
                                       value = Cw,
                                       method = rep("rep",lCw)))
    }
    return(events)
  })
}

#==================
#5. Custom function 
#==================
custom.func <- function(){
  return()
}

#==================
#6. ODEs System 
#==================
ode.func <- function(time, Initial.values, Parameters, custom.func){
  with( as.list(c(Initial.values, Parameters)),{
    
    dC <- ku*Cw - ke*C
    dCw <- 0 
    
    D.Magna <- C 
    
    list(c(dC=dC, dCw=dCw),"D.Magna" = D.Magna, "C_water" = Cw)
  })
}

#=================
# User input 
#=================
material <- "S2"  # Choose one from T1,T2,H1,S1,H2,S2
C_water <- c(1,0,10, 5) # NP concentration in water, in mg/L. Experimentally
# tested values: 0.1, 1.0 & 10.0 mg/L
sampling_time <- c(0,5,10, 25) #sampling time, in hours
user_input <- list("material"=material, "C_water"=C_water, "sampling_time" = sampling_time )

# Select a subset of the ODE output state variables to present on the UI
predicted.feats <- c("D.Magna",  "C_water")

# Create the parameters, initial values and events based on the user input
params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

simulation_time <- seq(0,30,1) # in hours
my_solution <-  deSolve::ode(times = simulation_time,  func = ode.func, y = inits, parms = params, 
                  custom.func = custom.func, method="lsodes",  events = events,rtol = 1e-6, atol = 1e-6)
print(head(my_solution[, c("time", "D.Magna",  "C_water")]))

###########################################################################################################################################
###  There are two options for connecting with the Jaqpot server:                                                                   #######
###   * Connect using the login.api() function which requests the API key, which can be downloaded from https://app.jaqpot.org/     #######
###   * Connect using the login.cred() function which requests the user's Jaqpot credentials (username & password)                  #######
## Please unhash the login.cred() line and hash the login.api() line if you want to login using the API key                         #######
###########################################################################################################################################

jaqpotr::login.cred()
#jaqpotr::login.api()

###########################################################################################################################################
###  Upload the biokinetics model on Jaqpot using the deploy.pbpk() function.                                                          ####
###  The function asks you (in the console) to name the model and provide a short description.                                         ####
###  Note that the model name cannot be changed, while the model description can be revised through the model's User Interface.        ####
###########################################################################################################################################

jaqpotr::deploy.pbpk(user_input, predicted.feats, create.params, create.inits, create.events, custom.func, ode.func) 

###########################################################################################################################################
###  Now let's try calling from here the model that lives in the Jaqpot server and acquire predictions, which we will then plot.        ###
###  The first step is to create a dataframe containing the input we want to provide to the model.                                      ###
###  If we don't remember the names of the independent features of the model, we can retrieve them from Jaqpot.                         ###
###  Note that you have to replace the modelID with the ID of your model which is printed above                                         ###
###########################################################################################################################################

# Set the model ID of the uploaded Jaqpot model
modelID <- "7KNh3EW3JoQEg33fE7vn"
model_feats <-  jaqpotr::get.model.feats(modelID = modelID) #function that retrieves independent feature names of the model
print(model_feats)

###########################################################################################################################################
###  Now let's create a dataframe with the names of the independent features and the corresponding values. After that, by providing     ###
###  the modelID and the dataframe we acquire predictions from our model that now "lives" in the Jaqpot server!                         ###
###########################################################################################################################################

# Create a dataframe containing the independent features of the model
# First we add all scalar model features to a dataframe
df =data.frame(sim.end = 30, sim.step = 1, sim.start = 0, material = "S2")
# The second step is to define the array features as lists on the already-created dataframe
df$C_water <-  list( c(1,0,10, 5))
df$sampling_time <- list(c(0,5,10, 20))
# Set the model ID of the uploaded Jaqpot model
modelID <- "7KNh3EW3JoQEg33fE7vn"
# Acquire model predictions
predictions <- jaqpotr::jaqpot.predict( df = df, modelID = modelID)

###########################################################################################################################################
###  Let's see what did the Jaqpot server return..                                                                                      ### 
###########################################################################################################################################

print(predictions)

###########################################################################################################################################
###  Finally, let's plot the D.Magna NM concentration alonside the NM concentration in water!                                           ###                                             ### 
###########################################################################################################################################

library(ggplot2)
plot_df <-  predictions$predictions
p1 = ggplot() + 
  geom_line(data = plot_df, aes(x = time, y = D.Magna), color = "blue") +
  xlab('time(hours)') +
  ylab('NM Concentration in D.magna (mg/g Wet Tissue)')

p2 = ggplot() + 
  geom_line(data = plot_df, aes(x = time, y = C_water), color = "red") +
  xlab('time(hours)') +
  ylab('NM Concentration (mg/L)')

print(p1)
print(p2)

