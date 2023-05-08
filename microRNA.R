### Set Directory ###
rm(list = ls())
library("data.table")

### Reading miRNA files and preprocess it ###
Protein_miRNA <- function(x) {
  z <- read.csv(x)
  y <- z$miRNAname
  y <- as.data.frame(y)
  colnames(y) <- "microRNAs"
  return(y)
}

miRNA <- function(x,y,z){
  S <- Protein_miRNA(y)
  M <- Protein_miRNA(z)
  x <- rbind(S,M)
  return(x)
}

### Calling Protein-miRNA function and saving the results ###  
setwd("~/Desktop/New")
x1 <- miRNA(x1, "HMOX1S.csv", "HMOX1M.csv")
x2 <- miRNA(x2, "INSRS.csv","INSRM.csv")
x3 <- Protein_miRNA("NOS3M.csv")
x4 <- miRNA(x4, "PPARGS.csv", "PPARGM.csv")
x5 <- miRNA(x5, "SLC2A4S.csv","SLC2A4M.csv")
x6 <- miRNA(x6, "SOD2S.csv", "SOD2M.csv")
x7 <- miRNA(x7, "TNFS.csv", "TNFM.csv")

setwd("~/Desktop/Project 2/MIRNA/P Proteins miRNAs/")
x8 <- Protein_miRNA("APOE.csv")
x9 <- Protein_miRNA("APOC1.csv")
x10 <- Protein_miRNA("HLA-C.csv")
x11 <- Protein_miRNA("HLA-DQA1.csv")
x12 <- Protein_miRNA("HLA-DRB1.csv")
x13 <- Protein_miRNA("NECTIN2.csv")
x14 <- Protein_miRNA("TOMM40.csv")

setwd("~/Desktop/Project 2/MIRNA")
write.csv(x1, file = "HMOX1F.csv", row.names = FALSE)
write.csv(x2, file = "INSRF.csv", row.names = FALSE)
write.csv(x3, file = "NOS3F.csv", row.names = FALSE)
write.csv(x4, file = "PPARGF.csv", row.names = FALSE)
write.csv(x5, file = "SLCA4F.csv", row.names = FALSE)
write.csv(x6, file = "SOD2F.csv", row.names = FALSE)
write.csv(x7, file = "TNFF.csv", row.names = FALSE)
write.csv(x8, file = "APOEF.csv", row.names = FALSE)
write.csv(x9, file = "APOC1F.csv", row.names = FALSE)
write.csv(x10, file = "HLA-CF.csv", row.names = FALSE)
write.csv(x11, file = "HLA-DQA1F.csv", row.names = FALSE)
write.csv(x12, file = "HLA-DRB1F.csv", row.names = FALSE)
write.csv(x13, file = "NECTIN2F.csv", row.names = FALSE)
write.csv(x14, file = "TOMM40F.csv", row.names = FALSE)
### Binding all miRNAs in one file ###
x89 <- rbind(x8,x9)
x1011 <- rbind(x10,x11)
x1213 <- rbind(x12,x13)
x891011 <- rbind(x89,x1011)
x121314 <- rbind(x1213,x14)
xf <- rbind(x891011,x121314)

Umirnas <- as.matrix(unique(xf))
flen <- length(Umirnas)


### Reading each miRNA files ###
miR <- function(x) {
  z <- read.csv(x)
  y <- as.matrix(unique(z$microRNAs))
  colnames(y) <- "microRNA"
  return(y)
}

setwd("~/Desktop/Project 2/MIRNA/Final miRNAs/")
HMOX1 <- miR("HMOX1F.csv")
lenHMOX1 <- length(HMOX1)
INSR <- miR("INSRF.csv")
lenINSR <- length(INSR)
NOS3 <- miR("NOS3F.csv")
lenNOS3 <- length(NOS3)
PPARG <- miR("PPARGF.csv")
lenPPARG <- length(PPARG)
SLC2A4 <- miR("SLCA4F.csv")
lenSLC2A4 <- length(SLC2A4)
SOD2 <- miR("SOD2F.csv")
lenSOD2 <- length(SOD2)
TNF <- miR("TNFF.csv")
lenTNF <- length(TNF)
APOE <- miR("APOEF.csv")
lenAPOE <- length(APOE)
APOC1 <- miR("APOC1F.csv")
lenAPOC1 <- length(APOC1)
TOMM40<- miR("TOMM40F.csv")
lenTOMM40 <- length(TOMM40)
NECTIN2<- miR("NECTIN2F.csv")
lenNECTIN2 <- length(NECTIN2)
HLA_C<- miR("HLA-CF.csv")
lenHLA_C <- length(HLA_C)
HLA_DQA1<- miR("HLA-DQA1F.csv")
lenHLA_DQA1 <- length(HLA_DQA1)
HLA_DRB1<- miR("HLA-DRB1F.csv")
lenHLA_DRB1 <- length(HLA_DRB1)

Adjacency_matrix <- matrix(data = 0, nrow = flen,
                           ncol = 8,byrow = FALSE)
colnames(Adjacency_matrix) <- c("microRNA", "APOE", 
                                "APOC1", "TOMM40","NECTIN2",
                                "HLA-C", "HLA-DQA1", "HLA-DRB1")

for (i in 1:flen) {
  Adjacency_matrix[i,1] = Umirnas[i,1]
}


fill <- function(t,y,w) {
  for (i in 1:flen) {
    m <- 0
    for (j in 1:lenHLA_DRB1) {
      if(Umirnas[i,1] == HLA_DRB1[j,1]){
        m <- 1
      }
      if(m == 1){
        Adjacency_matrix[i,8] <- 1
      }
    }  
  }
  return(Adjacency_matrix)
}

Adjacency_matrix[,2] <- fill(APOE,lenAPOE, 2)
Adjacency_matrix[,3] <- fill(APOC1,lenAPOC1, 3)
Adjacency_matrix[,4] <- fill(TOMM40,lenTOMM40, 4)
Adjacency_matrix[,5] <- fill(NECTIN2,lenNECTIN2, 5)
Adjacency_matrix[,6] <- fill(HLA_C,lenHLA_C, 6)
Adjacency_matrix[,7] <- fill(HLA_DQA1,lenHLA_DQA1, 7)
Adjacency_matrix[,8] <- fill(HLA_DRB1,lenHLA_DRB1, 8)


write.table(Adjacency_matrix, file = "adjacency_matrix.txt",row.names = FALSE,col.names =TRUE, sep = "\t" )
