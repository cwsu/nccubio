# read PAM1 from data
pam1<-read.table("C:/Users/User/��s��/�ͪ���T/pam1.txt")

# check PAM1 data
dim(pam1)
str(pam1)

# ��pam1�i��B�z
deal_pam1<-as.matrix(pam1)
deal_pam1<-deal_pam1/10000

# construct PAM250 from PAM1
recursive.matrixPower <- function(matrix,x) {
   if (x == 0)    return (matrix)
   else           return (matrix %*% recursive.matrixPower(matrix,x-1))
}
pam250<-recursive.matrixPower(deal_pam1,250)
pam250<-pam250*100

# output PAM250 as a file
write.table(pam250,file="C:/Users/User/��s��/�ͪ���T/pam250.txt")