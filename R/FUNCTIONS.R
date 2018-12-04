library(RcppAnnoy)
library(Matrix)
require(future)
require(future.apply)
library(rARPACK)
require(nmslibR)
require(data.table)
reticulate::use_python(system("dirname `which python3`",intern=T))

options(future.globals.maxSize= Inf)

#Diffusion Maps
#Based on Angerer, Philipp, et al. "destiny: diffusion maps for large-scale single-cell data in R." Bioinformatics 32.8 (2015): 1241-1243.
#With some modifications to make it fast.
ComputeDMap<-function(ex_mat,
                      nNN         = 10,
                      k           = 100,
                      nLocalsigma = 5L,
                      nThreads    = availableCores(),
                      nTrees      = 50,
                      M           = 10,
                      efC         = 100,
                      efS         = 100,
                      AnnMethod="Annoy",
                      EigDecompMethod="Irlba",
                      ...){
  plan(list(tweak(multicore, workers = nThreads)))
  
  f <- ncol(ex_mat)
  n <- nrow(ex_mat)
  
  start.time <- Sys.time()
  #knn
  if(AnnMethod=="Annoy"){
  
  assign("a", new(AnnoyEuclidean, f), envir = .GlobalEnv)
  print("Adding Items")
  for (i in seq(n)) {
  #  if(i %% 10000==0) {print(i)}
    a$addItem(i-1,ex_mat[i,])
  }
  print("Building Trees")
  a$build(nTrees)
  gc()
  
  NN.ind<-matrix(nrow=n,ncol=nNN)
  NN.dist<-matrix(nrow=n,ncol=nNN)
  ##################################################################
  #nNN<<-nNN
  fun.1<-function(x){
    a$getNNsByItem(x-1, nNN+1)[-1]
  }
  fun.2<-function(i){
    return(sapply(NN.ind[i,], function(x) a$getDistance(i-1,x)))
  }
  print("Getting NNs")
  NN.ind<-matrix(unlist(future_lapply(seq_len(n),fun.1),use.names=FALSE), ncol = nNN, byrow = TRUE)
  NN.dist<-matrix(unlist(future_lapply(seq_len(n),fun.2),use.names=FALSE), ncol = nNN, byrow = TRUE)
  } 
  if(AnnMethod=="Nmslib"){
    index_params = list('M'= M, 'indexThreadQty' = nThreads, 'efConstruction' = efC, 
                        'post' = 0, 'skip_optimized_index' = 1)
    query_time_params = list('efSearch' = efS)
    space_name = 'l2'     
    print("Building Index")
    init_nms <- NMSlib$new(input_data = ex_mat, Index_Params = index_params, 
                           Time_Params = query_time_params, space = space_name, 
                           space_params = NULL, method = 'hnsw', 
                           index_filepath = NULL, print_progress = TRUE)
    print("Getting NNs")
    tmp <- init_nms$knn_Query_Batch(ex_mat, k = nNN, num_threads = nThreads)
    NN.ind<-tmp$knn_idx-1
    NN.dist<-tmp$knn_dist
    rm(tmp)
  }
  gc()
  ###################################################################
  #sigma local 
  print("Computing local sigma")
  N=nLocalsigma
  sigma<-NN.dist[,N]/2
  trans.prob<-(NN.dist)^2
  
  vect_i<-rep(1:(n),nNN)
  vect_j<-as.vector(NN.ind)+1
  trans.prob<-as.vector(trans.prob)
  system.time(vect_x<-unlist(future_lapply(1:length(vect_i),
                                      function(x) sqrt(2*sigma[vect_i[x]]*sigma[vect_j[x]]/
                                                         (sigma[vect_i[x]]^2+sigma[vect_j[x]]^2))*
                                        exp(-trans.prob[x]/(sigma[vect_i[x]]^2+sigma[vect_j[x]]^2)))))
  tmp<-c(vect_i,vect_j)
  vect_j<-c(vect_j,vect_i)
  vect_i<-tmp
  rm(tmp)
  vect_x<-c(vect_x,vect_x)
  trans.prob<-sparseMatrix(i=vect_i-1,j=vect_j-1,x=vect_x,dims=c(n,n),index1=FALSE)
  rm(vect_x,vect_j,vect_i)
  
  gc()
  print("Computing transition probability matrix")
  diag(trans.prob)<-0
  trans.prob <- drop0(trans.prob)
  
  d <- rowSums(trans.prob)+1
  trans_p <- as(trans.prob, 'dgTMatrix')
  trans.prob<-sparseMatrix(trans_p@i, 
                           trans_p@j, 
                           x = trans_p@x / (d[trans_p@i + 1] * d[trans_p@j + 1]), 
                           dims = dim(trans_p), 
                           index1 = FALSE)
  
  rot<-Diagonal(x=(rowSums(trans.prob))^(-.5))
  
  gc()
  
  trans<-as(rot %*% trans.prob %*% rot, 'dgCMatrix')
  
  print("Computing eigendecomposition")
  k=k+1
  if(EigDecompMethod=="ARPACK"){decomp<-rARPACK::eigs(A=trans,k=k)}
  if(EigDecompMethod=="Irlba"){decomp<-irlba::partial_eigen(trans,n=k,tol=1e-5,work=100,symmetric = T)}
  end.time <- Sys.time()
  
  print(paste0("Total Computation Time: ",as.numeric(end.time-start.time, units = "secs")," secs"))
  
  return(list(vectors=as.matrix(t(t(decomp$vectors) %*% rot))[,-1],
              values =decomp$values))   
}

#Get approximate nearest neighbors

GetANNs<-function(mat,
                  nNN=10,
                  nThreads=availableCores(),
                  returnDistance=TRUE,
                  nTrees      = 50,
                  M           = 10,
                  efC         = 100,
                  efS         = 100,
                  AnnMethod="Annoy",...){
  
  plan(list(tweak(multiprocess, workers = nThreads)))
  
  f <- ncol(mat)
  if(AnnMethod=="Annoy"){
  assign("a", new(AnnoyEuclidean, f), envir = .GlobalEnv)
  n <- nrow(mat)                             
  
  for (i in seq(n)) {
#    if(i %% 100000==0) {print(i)}
    a$addItem(i-1,mat[i,])
  }
  
  a$build(nTrees)                           	
  
  fun.1<-function(x){
    a$getNNsByItem(x-1, nNN+1)[-1]
  }
  fun.2<-function(i){
    return(sapply(NN.ind[i,], function(x) a$getDistance(i-1,x)))
  }
  
  print("Getting NNs")
  NN.ind<-matrix(nrow=n,ncol=nNN)
  NN.ind<-matrix(unlist(future_lapply(seq_len(n),fun.1),use.names=FALSE), ncol = nNN, byrow = TRUE)
  if(returnDistance) {
    NN.dist<-matrix(nrow=n,ncol=nNN)
    NN.dist<-matrix(unlist(future_lapply(seq_len(n),fun.2),use.names=FALSE), ncol = nNN, byrow = TRUE)
  } 
  NN.ind<-NN.ind+1
  }
  if(AnnMethod=="Nmslib"){
    index_params = list('M'= M, 'indexThreadQty' = nThreads, 'efConstruction' = efC, 
                        'post' = 0, 'skip_optimized_index' = 1)
    query_time_params = list('efSearch' = efS)
    space_name = 'l2'     
    print("Building Index")
    init_nms <- NMSlib$new(input_data = mat, Index_Params = index_params, 
                           Time_Params = query_time_params, space = space_name, 
                           space_params = NULL, method = 'hnsw', 
                           index_filepath = NULL, print_progress = TRUE)
    print("Getting NNs")
    tmp <- init_nms$knn_Query_Batch(mat, k = nNN, num_threads = nThreads)
    NN.ind<-tmp$knn_idx
    if(returnDistance) {NN.dist<-tmp$knn_dist}
    rm(tmp)
  }
  if(returnDistance){
  return(list(NN.ind=NN.ind,
              NN.dist=NN.dist))
  } else {
    return(list(NN.ind=NN.ind)) 
  }
  
}

#Compute FLE
#3D Forceatlas2 with modifications to visualize huge graphs
#Forceatlas2 algorithm: Jacomy, Mathieu, et al. "ForceAtlas2, a continuous graph layout algorithm for handy network 
#visualization designed for the Gephi software." PloS one 9.6 (2014): e98679.
GetFLE<-function(inputFilePATH,
                 toolkitPATH,
                 memmory             = "1g",
                 layout              = "fa_3d",
                 nsteps              = 1000,
                 nThreads            = 1,
                 scalingRatio        = 1,
                 seed                = 1,
                 barnesHutTheta      = 1.2,
                 barnesHutUpdateIter = 1,
                 updateCenter        = F,
                 barnesHutSplits     = 1,
                 restart             = F
                 ){
  output=paste0(tools::file_path_sans_ext(inputFilePATH),"_FLE")
  if(restart==T){system(paste0("mv ",output,".distances.txt ",output,".distances_prev.txt"))}
  Command=paste0("java -Xmx",memmory," -Djava.awt.headless=true -cp ",
                 toolkitPATH,"force-atlas-3d.jar:",
                 toolkitPATH,"gephi-toolkit-0.9.2-all.jar:",
                 " GraphLayout ",
                 " --input ",inputFilePATH,
                 " --output ",output,
                 " --format txt",
                 if(restart){paste0(" --coords ",output,".txt")},
                 " --layout ",layout,
                 " --nsteps ",nsteps,
                 " --threads ",nThreads,
                 " --barnesHutTheta ",barnesHutTheta,
                 " --barnesHutUpdateIter ", barnesHutUpdateIter,
                 if(updateCenter){" --updateCenter true"} else {" --updateCenter false"},
                 " --updateCenter ",updateCenter,
                 " --scalingRatio ",scalingRatio,
                 " --barnesHutSplits ",barnesHutSplits,
                 " --seed ",seed
  )
  system(Command,wait=T)
  if(restart==T){
    df=rbind(read.table(paste0(output,".distances_prev.txt"),header=T),read.table(paste0(output,".distances.txt"),header=T))
    df$step<-0:(nrow(df)-1)
    write.table(df,file=paste0(output,".distances.txt"),quote = F,row.names = F,col.names = T)
  }
}
