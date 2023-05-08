setwd("~/Desktop/Drug lncRNA Project//")
library(phangorn)
library(data.table)
library(ape)
library(psych)
library(caret)
library(mlbench)
library(FunChisq)


Drugt <- read.csv("PaDEL  copy.csv",
                   row.names = 1, header = TRUE)
Drugt <- as.data.frame(Drugt)

s <- vector()
j <- 1
for (i in 1:1444) {
  if (Drugt[1,i] == Drugt[2,i] & 
      Drugt[2,i] == Drugt[3,i] & 
      Drugt[3,i] == Drugt[4,i] &
      Drugt[4,i] == Drugt[5,i] &
      Drugt[5,i] == Drugt[6,i] &
      Drugt[6,i] == Drugt[7,i] &
      Drugt[7,i] == Drugt[8,i] &
      Drugt[8,i] == Drugt[9,i] &
      Drugt[9,i] == Drugt[10,i] &
      Drugt[10,i] == Drugt[11,i] &
      Drugt[11,i] == Drugt[12,i] &
      Drugt[12,i] == Drugt[13,i] &
      Drugt[13,i] == Drugt[14,i] &
      Drugt[14,i] == Drugt[15,i] &
      Drugt[15,i] == Drugt[16,i] &
      Drugt[16,i] == Drugt[17,i] &
      Drugt[17,i] == Drugt[18,i] &
      Drugt[18,i] == Drugt[19,i] &
      Drugt[19,i] == Drugt[20,i] &
      Drugt[20,i] == Drugt[21,i] &
      Drugt[21,i] == Drugt[22,i] &
      Drugt[22,i] == Drugt[23,i] &
      Drugt[23,i] == Drugt[24,i] &
      Drugt[24,i] == Drugt[25,i] &
      Drugt[25,i] == Drugt[26,i] &
      Drugt[26,i] == Drugt[27,i] &
      Drugt[27,i] == Drugt[28,i] &
      Drugt[28,i] == Drugt[29,i] &
      Drugt[29,i] == Drugt[30,i] &
      Drugt[30,i] == Drugt[31,i] &
      Drugt[31,i] == Drugt[32,i] &
      Drugt[32,i] == Drugt[33,i] &
      Drugt[33,i] == Drugt[34,i] &
      Drugt[34,i] == Drugt[35,i] &
      Drugt[35,i] == Drugt[36,i] &
      Drugt[36,i] == Drugt[37,i] &
      Drugt[37,i] == Drugt[38,i] &
      Drugt[38,i] == Drugt[39,i] &
      Drugt[39,i] == Drugt[40,i]
  ){
    s[j] <- i
    j <- j+1
    
  }
  
}

Drugt <- Drugt[, -c(s)]

chisq.test(data_frame$treatment, data_frame$improvement, correct=FALSE)
Pear <- matrix(data = 0,nrow = 1151, ncol = 1151,
               byrow = FALSE)

colnames(Pear) <- colnames(m)
rownames(Pear) <- colnames(m)

for (i in 1:382) {
  for (j in 1:382) {
    x <- m[,i]
    y <- m[,j]
    Pear [i,j] <- cor(x, y, method = "pearson")
  }
}  

write.table(Pear, file = "m.txt",row.names = TRUE,col.names = TRUE, sep = "\t")






Drugt <- as.matrix(Drugt)

correlationMatrix <- cor(Drugt[,1:1151])
write.table(Drugt , file = "Drugcornew.txt",row.names = TRUE,col.names = TRUE, sep = "\t")

# summarize the correlation matrix
print(correlationMatrix)
correlationMatrix <- read.table("correlation.txt", sep = "\t",
                                row.names = 1, header = TRUE)
correlationMatrix <- as.matrix(correlationMatrix)

# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.6)

# print indexes of highly correlated attributes
print(highlyCorrelated)
Drugt <- Drugt[,-c(highlyCorrelated)]

