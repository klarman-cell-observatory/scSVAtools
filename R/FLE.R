source('/scSVAtools/R/FUNCTIONS.R');
FLE_InputParameters<-jsonlite::read_json("/home/ForceDirectedLayout_InputParameters.json",simplifyVector = T)
if(FLE_InputParameters$path==""){
  file_name<-"/home/NNG_output.csv"
} else {
  file_name<-FLE_InputParameters$path     
}

start.time <- Sys.time()
FLE<-GetFLE(inputFilePATH       = file_name,
            toolkitPATH         = "/gephi/",
            memory              = FLE_InputParameters$memory,
            nsteps              = FLE_InputParameters$nsteps,
            nThreads            = FLE_InputParameters$nThreads,
            scalingRatio        = FLE_InputParameters$scalingRatio,
            seed                = FLE_InputParameters$seed,
            barnesHutTheta      = FLE_InputParameters$barnesHutTheta,
            barnesHutUpdateIter = FLE_InputParameters$barnesHutUpdateIter,
            barnesHutSplits     = FLE_InputParameters$barnesHutSplits,
            updateCenter        = FLE_InputParameters$updateCenter,
            restart             = FLE_InputParameters$restart)

end.time <- Sys.time()
print(paste0("Total FLE Computation Time: ",as.numeric(end.time-start.time, units = "secs")," secs"))
