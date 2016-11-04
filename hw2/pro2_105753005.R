######################################
# the reference code of program2 
######################################

######################################
# initial
######################################
# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript pro2_105753005.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}
# commend1 : Rscript pro2_105753005.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta
# commend2 : Rscript pro2_105753005.R --input test.fasta --score pam100.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<-args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-args[i+1]
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

print("PARAMETERS")
print(paste("input file         :", i_f))
print(paste("output file        :", o_f))
print(paste("score file         :", s_f))
print(paste("gap open penalty   :", g_o))
print(paste("gap extend penalty :", g_e))

######################################
# main
######################################
# your code

# make the file become a object
inputFile <- file(i_f,open="r")
ifLine <- readLines(inputFile)
close(inputFile)

# read score table
scoreFile <- read.table(s_f)

#unlist : extract to vactor from a strsplit list
seqA <- unlist(strsplit(paste("*", ifLine[2], sep = ""),split = ""))
seqB <- unlist(strsplit(paste("*", ifLine[4], sep = ""),split = ""))

matrixRow <- length(seqA)
matrixCol <- length(seqB)

# create a matrix which value is NA 
seqMatrix <- matrix(NA, matrixRow, matrixCol, byrow = FALSE, dimnames = list(seqA,seqB))

# initial default score
seqMatrix[1,1] = 0
for(i in c(2:matrixRow)){
  seqMatrix[i,1] = as.numeric(g_o) + as.numeric(g_e)*(i-2)
}
for(j in c(2:matrixCol)){
  seqMatrix[1,j] = as.numeric(g_o) + as.numeric(g_e)*(j-2)
}

# calculate matrix value
gapOpen <- TRUE
for(i in c(2:matrixRow)){
  for(j in c(2:matrixCol)){
    slash <- seqMatrix[i-1,j-1]+scoreFile[seqA[i],seqB[j]]
    straight <- seqMatrix[i-1,j]+as.numeric(g_o)
    horizontal <- seqMatrix[i,j-1]+as.numeric(g_o)
    seqMatrix[i,j] <- max(slash,straight,horizontal)
  }
}

# trace path
rowId <- 2
colId <- 2
finalSeqA = seqA[2]
finalSeqB = seqB[2]

while(rowId < matrixRow || colId < matrixCol){
  tempMax <- max(seqMatrix[rowId+1,colId+1] , seqMatrix[rowId,colId+1] , seqMatrix[rowId+1,colId])
  if(tempMax == seqMatrix[rowId+1,colId+1]){
    finalSeqA <- paste(finalSeqA, seqA[rowId+1], seq ="")
    finalSeqB <- paste(finalSeqB, seqB[colId+1], seq ="")
    rowId <- rowId+1
    colId <- colId+1 
  }else if(tempMax == seqMatrix[rowId,colId+1]){
    finalSeqA <- paste(finalSeqA, "-", seq ="")
    finalSeqB <- paste(finalSeqB, seqB[colId+1], seq ="")
    colId <- colId+1 
  }else{
    finalSeqA <- paste(finalSeqA, seqA[rowId+1], seq ="")
    finalSeqB <- paste(finalSeqB, "-", seq ="")
    rowId <- rowId+1
  }
}

# write file
write(c(ifLine[1],finalSeqA,ifLine[3],finalSeqB),file = o_f)
print(finalSeqA)
print(finalSeqB)