# ----------------------------  INSTRUCTIONS ######
# -------------------------------------------------------- 
# First of all:
#   install isocir from .zip file
# Then:
#   install from the CRAN the following packages: combinat, circular, TSP
# 
# Change the folder where the data are in the computer:
  
  inputfolder <- "D:/"

# Change the output folder where you want the file with the N estimations of the global circular order to be written:

  outputfolder <- "D:/PCH/backward/10genes/fkh2"
#outputfolder <- "D:/PCH/6genes/ace2" #ejemplo de salida en el que se ejecuta con 6 genes 
# teniendo ace2 como último eliminado, al no poner la barra al final el nombre de los archivos
# comienza por ace2 y así se tienen juntos todos los archivos del paso de 6 genes en la misma carpeta.

# Change in the function eq.test the argument N to the number of simulations wanted.

library(isocir)

datosP <- read.csv2(paste(inputfolder,"11PombeCy.csv",sep=""),header=TRUE, row.names=1)
datosH <- read.csv2(paste(inputfolder,"11HumanCy.csv",sep=""),header=TRUE, row.names=1)


#### Drop the symbol: # in the line wanted to be executed.

## FORMER DATA (PERIODICITY):
#datosC <- read.table(paste(inputfolder,"11CerevisiaeI.txt",sep=""),header=TRUE, row.names=1)

## NEW ORTHOLOGS:
#datosC <- read.table(paste(inputfolder,"11Cerevisiae_neworthologsI.txt",sep=""),header=TRUE, row.names=1)
datosC <- read.csv2(paste(inputfolder,"11CerevisiaeCy.csv",sep=""),header=TRUE, row.names=1)

kappasH <- c(1.7211343, 2.3395968, 26.776229, 2.4580805) #cyclebase6:
kappasP <- c(1.6349211, 1.5415776, 1.6062273, 1.0778214, 9.0779582,1.400741, 0.1425247,  26.86056, 2.5155315,  3.458664)
kappasC <- c(1.8028731,  0.841727, 10.471244, 26.640786, 8.8043238, 1.7899144)


dataset <- rbind(as.matrix(datosP),as.matrix(datosC),as.matrix(datosH))

ws <- c(kappasP,kappasC,kappasH)
populations <- c(rep(1,10), rep(2,6), rep(3,4))
#populations <- c(rep(1,10), rep(2,6))

# -------------------------------------------------------- 
#### DROPPING GENES
# write in the vector called "genes" the number of the columns where the genes to be analized are.
# be careful to load again the wholedata


genes <- c(1:5,7:11) # The former final set of 7 genes for the 3 species OLD paralogs

#genes <- c(2,3,4,5,6,7) #initial set forward

dataset <- dataset[,genes]

set.seed(1)
date()
result.test <- eq.test(dataset, populations, ws=ws, method="TSP", control.method="alpha3", output=outputfolder, coef=6, N=2000)
date()
genes
colnames(datosP)[genes]
outputfolder3 <- outputfolder
result.test3 <- result.test
#result.test1
#outputfolder1

