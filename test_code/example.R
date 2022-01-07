library(DDSA)

# user input
data_dir = "../test_data/" # data directory
phase_range = c(-10,20) # phase range
wavelength_range = c(3500,8000) #wavelength range


# generate a supernova summary list
data_list = LoadData(DataDir = data_dir,PhaseRange = phase_range,WavelengthRange = wavelength_range)

#train mean
mean_result = Mean_SEDTraining(SNeList = data_list$SNeList ,Omega = data_list$Omega,isCV = FALSE)

# train principal component
R = 2  #The number of component
pc_result = PC_SEDTraining(R = R,SNeList =  data_list$SNeList,Theta0 = mean_result$Theta0,Omega = data_list$Omega,isCV = FALSE)


