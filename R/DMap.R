source('/scsvatools/FUNCTIONS.R');
start.time <- Sys.time()
DiffusionMap_InputParameters<-jsonlite::read_json("/home/DiffusionMap_InputParameters.json",simplifyVector = T)

if(summary(file(DiffusionMap_InputParameters$path))$class  ==  "gzfile"){
          ex_mat  <-  as.matrix(fread(cmd = paste0('gunzip -cq ',DiffusionMap_InputParameters$path),nThread=DiffusionMap_InputParameters$nThreads))
          } else {
          ex_mat  <-  as.matrix(fread(input = DiffusionMap_InputParameters$path,nThread=DiffusionMap_InputParameters$nThreads))
          }
dmap<-ComputeDMap(ex_mat          = ex_mat,
                  nNN             = DiffusionMap_InputParameters$nNN,
                  k               = DiffusionMap_InputParameters$k,
                  nLocalsigma     = DiffusionMap_InputParameters$nLocalsigma,
                  nThreads        = DiffusionMap_InputParameters$nThreads,
                  nTrees          = DiffusionMap_InputParameters$nTrees,
                  M               = DiffusionMap_InputParameters$M,
                  efC             = DiffusionMap_InputParameters$efC,
                  efS             = DiffusionMap_InputParameters$efS,
                  AnnMethod       = DiffusionMap_InputParameters$AnnMethod,
                  EigDecompMethod = DiffusionMap_InputParameters$EigDecompMethod
                )
end.time <- Sys.time()
print(paste0("Total DMap Computation Time: ",as.numeric(end.time-start.time, units = "secs")," secs"))
fwrite(x         = as.data.frame(dmap$vectors),
       file      = "/home/DMap_eigenvectors_output.csv",
       nThread   = DiffusionMap_InputParameters$nThreads,
       quote     = F,
       row.names = F,
       col.names = F       
       )

fwrite(x         = as.data.frame(dmap$values),
       file      = "/home/DMap_eigenvalues_output.csv",
       quote     = F,
       row.names = F,
       col.names = F       
       )
