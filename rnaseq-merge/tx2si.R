#R script to calculate Splicing Index for Tx/Exon etc.


####
#Prerequisites
####

library(argparser,quietly =T)

####
#Version
####

version="0.3"

#v0.1, for tx SI calculation
#v0.2, add filter to reduce calculation. Change gene column to 3rd col
#v0.3, use parallel to speed up the calculation

description=paste0("tx2si\nversion ",version,"\n","Usage:\nDescription: Calculate Splicing Index for tx/exon\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)

parser <- add_argument(parser, arg="--gene",short="-g", type="character", help = "Gene Expr")
parser <- add_argument(parser, arg="--tx",short="-t", type="character", help = "Tx/Exon Expr")
parser <- add_argument(parser, arg="--anno",short="-a", type="character", help = "Annotation file for tx/exon")
parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Output, Splicing Index")
parser <- add_argument(parser, arg="--filter",type="character", help = "Count filter for all the exps",default="auto")

args = parse_args(parser)

print(args)

#other parameters


###########
#Libraries
##########

library(parallel)

###########
#Functions
###########

write_table_proper<-function(file,data,name="Gene") {
	data.df<-data.frame(name=rownames(data),data)
	names(data.df)<-c(name,colnames(data))

	write.table(data.df,file=file, row.names=FALSE,sep="\t",quote=F)

	#write.table(data.frame(name=rownames(data),data),file=file, row.names=FALSE,sep="\t",quote=F)
}

filter_data <- function(mat,type="sum",cutoff=10,na.rm=0) {
	#remove NA
	mat[is.na(mat)]<-na.rm


	if(cutoff=="auto") {
	  num.cutoff<-ncol(mat)*5 # count cutoff eq two times number of samples #changed to *5 3/26
	} else {
		if(grepl("^auto",cutoff,perl=T)) {
			fold<-regmatches(cutoff,regexpr("\\d+$",cutoff,perl=T))
			num.cutoff<-ncol(mat)*as.numeric(fold)
		} else {
			num.cutoff<-as.numeric(cutoff)
		}
	}

		
	if(type == "sum") {
		mat.sel<-mat[apply(mat,1,sum)>=num.cutoff,]
	} else if (type=="none") {
		mat.sel<-mat
	}
	
	return(mat.sel)
}


cal_si<-function(txname) {
	tx_expr_vec<-tx_expr.sel[txname,]
	genename<-anno[txname,2] #change gene name into 3rd column of anno
	if(genename %in% rownames(gene_expr)) {
		gene_expr_vec<-gene_expr[genename,]
		as.vector(unlist(tx_expr_vec))/as.vector(unlist(gene_expr_vec))
	} else {
		rep(0,length(tx_expr_vec))
	}
}


###########
#Program starts
##########

#if gene exp is 0, return inf?
#tx_addup<-0.001

gene_expr<-read.table(args$"gene",header=T,row.names=1,sep="\t",check.names=F,flush=T,comment.char="",quote="")
tx_expr<-read.table(args$"tx",header=T,row.names=1,sep="\t",check.names=F,flush=T,comment.char="",quote="")

cat(nrow(tx_expr)," features read for SI calculation.\n")

#anno first column as tx, second col as gene
anno<-read.table(args$anno,header=T,row.names=1,sep="\t",check.names=F,flush=T,comment.char="",quote="")

#save image
rdatafile=sub(".txt$",".rdata",args$out,perl=T)


#filter data to reduce calculation
tx_expr.sel<-filter_data(tx_expr,cutoff=args$filter)

cat(nrow(tx_expr.sel)," features used for SI calculation.\n")

#detect cores
ncores<-detectCores()
cat(ncores," cores detected\n")

#open a cluster
cl <- makeCluster(detectCores())
clusterExport(cl,list("tx_expr.sel","anno","gene_expr"))

tx_si <- t(parSapply(cl, rownames(tx_expr.sel),cal_si))

stopCluster(cl)


rownames(tx_si)<-rownames(tx_expr.sel)
colnames(tx_si)<-colnames(tx_expr.sel)

#
write_table_proper(tx_si,file=args$out,"Feature")

