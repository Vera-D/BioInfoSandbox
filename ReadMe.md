## Project Description - this file contains R script code for assesing sequence complexity with common genes used in oncology dx assays when looking at 1 mer,2,3,4,5,and 6 mers. The goal is to get an idea of sequence distriution for the different genes across the differenet n-mer combos. In this script we will have a FASTA file with the gene sequence of a panel. A simple panel will sequence 6 genes: KRAS, NRAS, HRAS. In the second run I will include BRAF, EGFR, PI3CK. To Test: Copy paste the code per section onto the R IDE. Be sure to modify the paths to your box. 

## Step1: Go to Ensemble and download the reference sequence per gene using CHR38 annotatations. 

## HRAS is small so lets start there

library("seqinr", lib.loc="~/R/win-library/3.3")

HRAS <-read.fasta(file="C:/Users/vdiaz/Desktop/BioInfo-Sandbox/Homo_sapiens_HRAS_sequence.fa", as.string = FALSE, seqtype = "DNA")
HRAS<-HRAS[[1]]
length(HRAS)
attr(HRAS,"Annot")
par(mfrow=c(3,2))

barplot(count(HRAS,1),col="forest green",main="HRAS count per base")
barplot(count(HRAS,2),col="forest green",main="HRAS count per 2-mer",las=2)
barplot(count(HRAS,3),col="forest green",main="HRAS count per 3-mer",las=2)
barplot(count(HRAS,4),col="forest green",main="HRAS count per 4-mer",las=2)
barplot(count(HRAS,5),col="forest green",main="HRAS count per 5-mer",las=2)
barplot(count(HRAS,6),col="forest green",main="HRAS count per 6-mer",las=2)

## Lets plot KRAS
KRAS <-read.fasta(file="C:/Users/vdiaz/Desktop/BioInfo-Sandbox/Homo_sapiens_KRAS_sequence.fa", as.string = FALSE, seqtype = "DNA")
KRAS<-KRAS[[1]]
length(KRAS)
attr(KRAS,"Annot")
par(mfrow=c(3,2))

barplot(count(KRAS,1),col="forest green",main="KRAS count per base")
barplot(count(KRAS,2),col="forest green",main="KRAS count per 2-mer",las=2)
barplot(count(KRAS,3),col="forest green",main="KRAS count per 3-mer",las=2)
barplot(count(KRAS,4),col="forest green",main="KRAS count per 4-mer",las=2)
barplot(count(KRAS,5),col="forest green",main="KRAS count per 5-mer",las=2)
barplot(count(KRAS,6),col="forest green",main="KRAS count per 6-mer",las=2)


## Lets plot TP53
TP53 <-read.fasta(file="C:/Users/vdiaz/Desktop/BioInfo-Sandbox/Homo_sapiens_TP53_sequence.fa", as.string = FALSE, seqtype = "DNA")
TP53<-TP53[[1]]
length(TP53)
attr(TP53,"Annot")
par(mfrow=c(3,2))

barplot(count(TP53,1),col="forest green",main="TP53 count per base")
barplot(count(TP53,2),col="forest green",main="TP53 count per 2-mer",las=2)
barplot(count(TP53,3),col="forest green",main="TP53 count per 3-mer",las=2)
barplot(count(TP53,4),col="forest green",main="TP53 count per 4-mer",las=2)
barplot(count(TP53,5),col="forest green",main="TP53 count per 5-mer",las=2)
barplot(count(TP53,6),col="forest green",main="TP53 count per 6-mer",las=2)

## Lets output some values
## Set the working directory
setwd("C:/Users/vdiaz/Desktop/BioInfo-Sandbox/KRASout")
write.csv(count(KRAS,1), file = "1mer.csv")
write.csv(count(KRAS,2), file = "2mer.csv")
write.csv(count(KRAS,3), file = "3mer.csv")
write.csv(count(KRAS,4), file = "4mer.csv")
write.csv(count(KRAS,5), file = "5mer.csv")
write.csv(count(KRAS,6), file = "6mer.csv")

setwd("C:/Users/vdiaz/Desktop/BioInfo-Sandbox/HRASout")
write.csv(count(HRAS,1), file = "1mer.csv")
write.csv(count(HRAS,2), file = "2mer.csv")
write.csv(count(HRAS,3), file = "3mer.csv")
write.csv(count(HRAS,4), file = "4mer.csv")
write.csv(count(HRAS,5), file = "5mer.csv")
write.csv(count(HRAS,6), file = "6mer.csv")

setwd("C:/Users/vdiaz/Desktop/BioInfo-Sandbox/TP53out")
write.csv(count(TP53,1), file = "1mer.csv")
write.csv(count(TP53,2), file = "2mer.csv")
write.csv(count(TP53,3), file = "3mer.csv")
write.csv(count(TP53,4), file = "4mer.csv")
write.csv(count(TP53,5), file = "5mer.csv")
write.csv(count(TP53,6), file = "6mer.csv")

## Probabilities
HRAS_4mer<-count(HRAS,4)
KRAS_4mer<-count(KRAS,4)
TP53_4mer<-count(TP53,4)
KRAS_6mer<-count(KRAS,6)
TP53_6mer<-count(TP53,6)
HRAS_6mer<-count(HRAS,6)

## PROBABILIY CALC
P_TP53_4mer<-TP53_4mer/256
P_HRAS_4mer<-HRAS_4mer/256
P_KRAS_4mer<-KRAS_4mer/256
P_TP53_6mer<-TP53_6mer/4096
P_KRAS_6mer<-KRAS_6mer/4096
P_HRAS_6mer<-HRAS_6mer/4096

## as fractions
fractions(P_TP53_4mer)

## Function to count how many times get 0 in the vector
func1 <- function(vector){
counter=0
	for (i in 1:length(vector)) { 
		if(vector[i] == 0) {	
			counter = counter+1
		}
    }
	return(counter)
}

