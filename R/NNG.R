source('/scSVAtools/R/FUNCTIONS.R');

NNG_InputParameters<-jsonlite::read_json("/home/NNG_InputParameters.json",simplifyVector = T)

if(NNG_InputParameters$path==""){
  file_name<-"/home/DMap_eigenvectors_output.csv"
} else {
  file_name<-NNG_InputParameters$path     
}

start.time <- Sys.time()
mat<-as.matrix(fread(file_name,nThread = NNG_InputParameters$nThreads))
NNs<-GetANNs(mat,
             nNN             = NNG_InputParameters$nNN,
             nThreads        = NNG_InputParameters$nThreads,
             returnDistance  = FALSE,
             nTrees          = NNG_InputParameters$nTrees,
             M               = NNG_InputParameters$M,
             efC             = NNG_InputParameters$efC,
             efS             = NNG_InputParameters$efS,
             AnnMethod       = NNG_InputParameters$AnnMethod
            )

end.time <- Sys.time()
print(paste0("Total NNG Computation Time: ",as.numeric(end.time-start.time, units = "secs")," secs"))
#save adjacency table

fwrite(x         = as.data.frame(cbind(1:nrow(NNs$NN.ind),NNs$NN.ind)),
       file      = "/home/NNG_output.csv",
       nThread   =  NNG_InputParameters$nThreads,
       quote     = F,
       row.names = F,
       col.names = F
       )
