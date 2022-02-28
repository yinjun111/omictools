#R script to calculate Splicing Index for Tx


####
#Prerequisites
####

library(argparser,quietly =T)

####
#Version
####

version="0.1"

description=paste0("tx2si\nversion ",version,"\n","Usage:\nDescription: Calculate Splicing Index for tx\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)

parser <- add_argument(parser, arg="--gene",short="-g", type="character", help = "Gene TPM")
parser <- add_argument(parser, arg="--tx",short="-t", type="character", help = "Tx TPM")
parser <- add_argument(parser, arg="--anno",short="-a", type="character", help = "Annotation file for tx")
parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Output, Splicing Index")

args = parse_args(parser)

print(args)

#other parameters


###########
#Libraries
##########


###########
#Functions
###########

write_table_proper<-function(file,data,name="Gene") {
	data.df<-data.frame(name=rownames(data),data)
	names(data.df)<-c(name,colnames(data))

	write.table(data.df,file=file, row.names=FALSE,sep="\t",quote=F)

	#write.table(data.frame(name=rownames(data),data),file=file, row.names=FALSE,sep="\t",quote=F)
}


###########
#Program starts
##########

tx_addup<-0.001

gene_tpm<-read.table(args$"gene",header=T,row.names=1,sep="\t",check.names=F,flush=T,comment.char="",quote="")
tx_tpm<-read.table(args$"tx",header=T,row.names=1,sep="\t",check.names=F,flush=T,comment.char="",quote="")

#anno first column as tx, second col as gene
anno<-read.table(args$anno,header=T,row.names=1,sep="\t",check.names=F,flush=T,comment.char="",quote="")


#save image
rdatafile=sub(".txt$",".rdata",args$out,perl=T)


tx_si<-c()

txs<-c()

#calculate SI
for(rnum in 1:nrow(tx_tpm)) {
	genename<-anno[rownames(tx_tpm)[rnum],1]
	if(genename %in% rownames(gene_tpm)) {
		tx_si<-rbind(tx_si,tx_tpm[rnum,]/gene_tpm[genename,])
		txs<-c(txs,rownames(tx_tpm)[rnum])
	}
}

rownames(tx_si)<-txs
colnames(tx_si)<-colnames(tx_tpm)

#
write_table_proper(tx_si,file=args$out,"Feature")







