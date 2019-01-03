require(rhdf5)

#' Generate hdf5 input file from Seurat object
#'
#' @param SeuratObject Name of the Seurat Object to be saved in hdf5 format   
#' @param FileName     Output file name (with extension .h5,hdf5,loom)
#' @param ChunkSize    Size of the chunks
#' @param CompressionLevel Compression Level
#' @return Seurat Object saved in hdf5 file format
#' @examples 
#' GenerateInputFiles_fromSeurat(PBMC,"PBMC.h5",ChunkSize=1000)
#' @export

GenerateInputFiles_fromSeurat<-function(SeuratObject,FileName,ChunkSize=1000,CompressionLevel=9){
  h5createFile(FileName)  
  DimExMat<-dim(SeuratObject@data)
  h5createDataset(file=FileName, "matrix", c(DimExMat[1],DimExMat[2]),
                  storage.mode = "double", chunk=c(ChunkSize,ChunkSize), level=CompressionLevel)
  print("Saving Gene Expression Matrix")
  d<-1:DimExMat[2];n=ChunkSize;chunks<-split(d, ceiling(seq_along(d)/n))
  for(i in 1:length(chunks)){
    print(i)  
    mat<-as.matrix(SeuratObject@data[,chunks[[i]]])
    h5write(mat, 
            file=FileName,
            name="matrix", 
            index=list(NULL,min(chunks[[i]]):max(chunks[[i]])))
  }
  #GeneNames
  print("Saving Row Attributes")
  h5createGroup(FileName,"row_attrs")
  h5createDataset(file=FileName, 
                  "row_attrs/GeneNames",
                  dims = nrow(SeuratObject@data),
                  size = max(nchar(rownames(SeuratObject@data)))+1L,
                  chunk = ChunkSize,
                  storage.mode = "character", 
                  level=9)
  h5write(rownames(SeuratObject@data), 
          file=FileName,
          name="row_attrs/GeneNames")
  #MetaData
  print("Saving Column Attributes")
  Group<-"col_attrs"
  h5createGroup(FileName,Group)
  meta.data<-SeuratObject@meta.data
  emb<-names(SeuratObject@dr)
  
  for(i in emb){
  meta.data<-cbind(meta.data,SeuratObject@dr[[i]]@cell.embeddings)
  }

  meta.data.class<-sapply(meta.data, class)
  
  Ncol<-nrow(meta.data)
  for(i in 1:length(meta.data.class)){
    
    DataSetName <- paste0(Group,"/",names(meta.data.class)[i])
    DataSetType <- meta.data.class[i]
    
    if(DataSetType=='numeric') {
      
      h5createDataset(file=FileName, DataSetName,dims = Ncol,
                      storage.mode = "double", level=9)
      h5write(meta.data[,i], file=FileName,name=DataSetName)
      
    }
    if(DataSetType=='integer') {
      h5createDataset(file=FileName, DataSetName,dims = Ncol,
                      storage.mode = "integer", level=9)
      h5write(meta.data[,i], file=FileName,name=DataSetName)
    }
    if(DataSetType=='character') {
      h5createDataset(file=FileName, DataSetName,dims = Ncol,
                      storage.mode = "integer", level=9)
      tmp=as.factor(meta.data[,i])
      h5write(as.integer(tmp), file=FileName,name=DataSetName)
      
      l<-levels(tmp)
      h5createDataset(file=FileName, 
                      paste0(DataSetName,"_l"),
                      dims = length(l),
                      size = length(l),
                      storage.mode = "character", 
                      level=9)
      h5write(l, 
              file=FileName,
              name=paste0(DataSetName,"_l"))
    }
    if(DataSetType=='factor') {
      h5createDataset(file=FileName, DataSetName,dims = Ncol,
                      storage.mode = "integer", level=9)
      h5write(as.integer(meta.data[,i]), file=FileName,name=DataSetName)
      
      l<-levels(meta.data[,i])
      h5createDataset(file=FileName, 
                      paste0(DataSetName,"_l"),
                      dims = length(l),
                      size = max(nchar(l))+1L,
                      storage.mode = "character", 
                      level=9)
      h5write(l, 
              file=FileName,
              name=paste0(DataSetName,"_l"))
    }
  }
  H5close()
}

