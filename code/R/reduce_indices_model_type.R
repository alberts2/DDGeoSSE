for (i in 1:400000){
  test = read.csv(paste0("sim.",i,".labels.csv"))
  test$model_type = test$model_type-1
  write.csv(test, paste0("../simulate2/sim.",i,".labels.csv"))
}