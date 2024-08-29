set.seed(1234)
#library(dplyr)
#library(ggplot2)
library(tidyverse)
library(curl)
#library(wordspace) # for nromalizing rows or columns in unit length
#library(Matrix)
library(DAAG)
library(keras)
library(tensorflow)
library(doParallel)
intdata_big_normalize<-read.table("/scratch/iparvez/analysis_single/intdata_big_normalize.txt")
data<-t(intdata_big_normalize)

#Assigining number of cores
cl <- makeCluster(20)
registerDoParallel(cl)
tf$random$set_seed(1234)

# set model
model <- keras_model_sequential()
model %>%
  layer_dense(units = 500, activation = "tanh", kernel_regularizer = regularizer_l1_l2(l1 = 0.0000001, l2 =0.9 ),input_shape = ncol(data),kernel_initializer=initializer_random_uniform(minval = -0.0001, maxval = 0.0001,seed = 1234)) %>%
  layer_dropout(rate=0.5,seed=1234) %>%
  layer_dense(units = 150, activation = "tanh", kernel_regularizer = regularizer_l1_l2(l1 = 0.0000001, l2 =0.9 ),name = "bottleneck",kernel_initializer=initializer_random_uniform(minval = -0.0001, maxval = 0.0001,seed = 1234)) %>%
  layer_dropout(rate=0.5,seed=1234) %>%
  layer_dense(units = 500, activation = "tanh",kernel_regularizer = regularizer_l1_l2(l1 = 0.0000001, l2 =0.9 ),kernel_initializer=initializer_random_uniform(minval = -0.0001, maxval = 0.0001,seed = 1234)) %>%
  layer_dropout(rate=0.5,seed=1234) %>%
  layer_dense(units = ncol(data),kernel_regularizer = regularizer_l1_l2(l1 = 0.0000001, l2 =0.9 ),kernel_initializer=initializer_random_uniform(minval = -0.0001, maxval = 0.0001,seed = 1234))



# view model layers
summary(model)

# compile model
model %>% compile(
  loss = "mean_squared_error",
  optimizer = optimizer_adam(learning_rate=0.0001)
)




# fit model
model %>% fit(
  x = data,
  y = data,
  epochs = 20,
  use_multiprocessing=TRUE,
  batch_size = dim(data)[1]/10,
  verbose = 1
)



# evaluate the performance of the model
mse.ae2 <- evaluate(model, data, data)
# mse.ae2



stopCluster(cl)
# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
auto_500150500_20_all<- predict(intermediate_layer_model, data)
write.table(auto_500150500_20_all,"/scratch/iparvez/analysis/auto_500150500_20_all.txt")
write.table(mse.ae2,"/scratch/iparvez/analysis/mse_auto_500150500_20_all.txt")
