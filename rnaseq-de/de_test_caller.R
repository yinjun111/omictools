

####
#Prerequisites
####

library(argparser,quietly =T)

####
#Version
####

version="0.8"

#0.2b, change auto filter to *5. Add indfilter and cookscutoff option
#0.23, add write_table_proper
#0.3, add MA plot,fixed xlim for volcano plot
#0.31, add plotting option
#0.4 add new volcano plot and gene annotation. --anno and --config are used. --anno is different from previous argument
#0.41, changed title for Significance 
#0.51, add xlim ylim for volcano
#0.6, supports complicated GLM analysis. Change formula into model.matrix
#0.61, minor changes for read.table
#0.62, remove model matrix empty columns
#0.63, change library calling lcoations
#0.64, fix bug for group names starting with numbers
#0.7, add LM based analysis for Splicing Index etc.
#0.71, change significance column name
#0.72, add inf removal
#0.8, add NOISeq to support no replicate

description=paste0("de_test\nversion ",version,"\n","Usage:\nDescription: Differential Expression calculation\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)

parser <- add_argument(parser, arg="--in",short="-i", type="character", help = "Expr input file")
parser <- add_argument(parser, arg="--config",short="-c", type="character", help = "Config file for samples")
parser <- add_argument(parser, arg="--anno",short="-a", type="character", help = "Annotation file for genes")
parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Output file")
parser <- add_argument(parser, arg="--formula",short="-f", type="character", help = "DESeq formula")
parser <- add_argument(parser, arg="--treat",short="-t", type="character", help = "treatment name")
parser <- add_argument(parser, arg="--ref",short="-r", type="character", help = "reference name")
parser <- add_argument(parser, arg="--fccutoff", type="float", help = "Log2 FC cutoff",default=1)
parser <- add_argument(parser, arg="--qcutoff", type="float", help = "qcutoff",default=0.05)
parser <- add_argument(parser, arg="--pmethod",type="character", help = "Method used to calculate P value",default="DESeq2-Wald")
parser <- add_argument(parser, arg="--qmethod",type="character", help = "FDR Method",default="BH")
parser <- add_argument(parser, arg="--filter",type="float", help = "Count filter for all the exps. Default as mininum of (No. of samples) x 5 counts",default="auto")
parser <- add_argument(parser, arg="--independentfiltering",type="logical", help = "DESeq2 independentFiltering",default=T)
parser <- add_argument(parser, arg="--cookscutoff",type="logical", help = "DESeq2 cooksCutoff",default=T)
parser <- add_argument(parser, arg="--plot",type="logical", help = "Whether to generate volcano and MA plots",default=T)

args = parse_args(parser)

print(args)

#other parameters

core=1

#args_file = "tempArgObjectFile.rds"
#saveRDS(args, args_file); print(args); quit(); #comment this after creating args_file
#args = readRDS(args_file)  

#independentfiltering=T
#cookscutoff=T

###########
#Libraries
##########


library(DESeq2,quietly =T)
library(ggplot2,quietly =T)
library(NOISeq,quietly =T)
library(EnhancedVolcano,quietly =T)
library(Cairo,quietly =T)



#####
#Filter genes
#####


filter_data <- function(mat,type="sum",cutoff=10,na.rm=0) {
	#remove NA and Inf
	is.na(mat) <- sapply(mat, is.infinite)
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

	#changed here to return T/F	
	if(type == "sum") {
		#mat.sel<-mat[apply(mat,1,sum)>=num.cutoff,]
		return(apply(mat,1,sum)>=num.cutoff)
	} else if (type=="none") {
		#mat.sel<-mat
		return(rep(TRUE,nrow(mat)))
	}
	
	#return(mat.sel)
}


######
#DESeq2
######
deseq2_test <- function(mat,anno,design,fc_cutoff=1,q_cutoff=0.05,pmethod="Wald",qmethod="BH",core=core,treat,ref,independentfiltering=T,cookscutoff=T){

	#Get the last variable
	design<-gsub(" ","",design)
	design.vars<-gsub("~","",design)
	design.vars<-unlist(strsplit(design.vars,"\\+"))
	
	#change treat and ref names
	treat.m=make.names(treat)
	ref.m=make.names(ref)
	
	#anno relevel column for last variable #only works for 1vs1 for the selected col
	#convert to factor
	anno[,design.vars[length(design.vars)]]<-factor(make.names(anno[,design.vars[length(design.vars)]]),levels=c(ref.m,treat.m))
	
	#convert to model matrix
	design.mm<-model.matrix(as.formula(design),anno)	
	
	#remove empty col
	design.mm<-design.mm[,apply(design.mm,2,sum)>0] 
	print(design.mm)
	
	dds <- DESeqDataSetFromMatrix(countData = round(mat),
	                              colData = anno,
	                              design= design.mm)
	contrast<-paste0(design.vars[length(design.vars)],treat.m)
	
	#perform test
	dds <- DESeq(dds,test=pmethod)
	resultsNames(dds) # lists the coefficients
	
	#comp<-tail(all.vars(as.formula(design)),n=1) #the last variable of the formula is used as comparison
	
	if(cookscutoff) {
		res <- results(dds, pAdjustMethod =qmethod,contrast=list(contrast),independentFiltering=independentfiltering)
	} else {
		res <- results(dds, pAdjustMethod =qmethod,contrast=list(contrast),independentFiltering=independentfiltering,cooksCutoff=cookscutoff)	
	}
	
	fc<-res[,2]
	p<-res[,5]
	q<-res[,6]
	stat<-apply(res[,c(1,3,4)],1,function(x) {paste(x,collapse = ",")})
	
	#significance by fc & qval
	sig<-rep(0,length(q))
	
	sig[!is.na(fc) & fc>=fc_cutoff & !is.na(q) & q<q_cutoff]=1
	sig[!is.na(fc) & fc<=-fc_cutoff & !is.na(q) & q<q_cutoff]=-1
	sig[is.na(fc) | is.na(q) ]=NA

	mat.result<-cbind(fc,stat,p,q,sig)
	
	#colnames(mat.result)<-c(values(res)[[2]][2],"DESeq2 Stat:Mean,SE,Wald stat",values(res)[[2]][5],values(res)[[2]][6],paste("Significance: Abs(Log2FC)>=",round(fc_cutoff,3)," ",qmethod, "P<",q_cutoff,sep=""))
	
	colnames(mat.result)<-c(paste0("Log2FC ",treat," vs ",ref),"DESeq2 Stat:Mean,SE,Wald stat",values(res)[[2]][5],values(res)[[2]][6],paste("Significance: Abs(Log2FC)>=",round(fc_cutoff,3)," ",qmethod, "P<",q_cutoff,sep=""))
	
	result<-list()
	result$result=mat.result
	result$dds=dds
	result$res=res

	return(result)
}

deseq2_test_v1 <- function(mat,anno,design,fc_cutoff=1,q_cutoff=0.05,pmethod="Wald",qmethod="BH",core=core,treat,ref,independentfiltering=T,cookscutoff=T){

  #anno and design may need to be checked


	dds <- DESeqDataSetFromMatrix(countData = round(mat),
	                              colData = anno,
	                              design= as.formula(design))
	dds <- DESeq(dds,test=pmethod)
	resultsNames(dds) # lists the coefficients
	
	comp<-tail(all.vars(as.formula(design)),n=1) #the last variable of the formula is used as comparison
	
	if(cookscutoff) {
		res <- results(dds, pAdjustMethod =qmethod,contrast=c(comp,treat,ref),independentFiltering=independentfiltering)
	} else {
		res <- results(dds, pAdjustMethod =qmethod,contrast=c(comp,treat,ref),independentFiltering=independentfiltering,cooksCutoff=cookscutoff)	
	}
	
	fc<-res[,2]
	p<-res[,5]
	q<-res[,6]
	stat<-apply(res[,c(1,3,4)],1,function(x) {paste(x,collapse = ",")})
	
	#significance by fc & qval
	sig<-rep(0,length(q))
	
	sig[!is.na(fc) & fc>=fc_cutoff & !is.na(q) & q<q_cutoff]=1
	sig[!is.na(fc) & fc<=-fc_cutoff & !is.na(q) & q<q_cutoff]=-1
	sig[is.na(fc) | is.na(q) ]=NA

	mat.result<-cbind(fc,stat,p,q,sig)
	
	colnames(mat.result)<-c(values(res)[[2]][2],"DESeq2 Stat:Mean,SE,Wald stat",values(res)[[2]][5],values(res)[[2]][6],paste("Significance: ",treat," vs ",ref," Abs(Log2FC)>=",round(fc_cutoff,3)," ",qmethod, "P<",q_cutoff,sep=""))
	
	
	
	result<-list()
	result$result=mat.result
	result$dds=dds
	result$res=res

	return(result)
}

######
#NOISeq for single rep
######

noiseq_test_norep <- function(mat,anno,design,fc_cutoff=1,q_cutoff=0.05,treat,ref){

	#noiseq for single replicate

	#Get the last variable
	design<-gsub(" ","",design)
	design.vars<-gsub("~","",design)
	design.vars<-unlist(strsplit(design.vars,"\\+"))
	
	#change treat and ref names
	treat.m=make.names(treat)
	ref.m=make.names(ref)
	
	#anno relevel column for last variable #only works for 1vs1 for the selected col
	#convert to factor
	anno[,design.vars[length(design.vars)]]<-factor(make.names(anno[,design.vars[length(design.vars)]]),levels=c(treat.m,ref.m))
	
	#create noiseq data
	mydata<-readData(data = mat,factors = anno)

	#de test
	myresults<-noiseq(mydata, factor = design.vars[length(design.vars)], k = NULL, norm = "tmm", replicates = "no")
	
	res<-myresults@results[[1]]
	
	fc<-res[,3]
	p<-1-res[,5]
	q<-p #in noiseq, p is q
	stat<-apply(res[,c(3,4,6)],1,function(x) {paste(x,collapse = ",")})
	
	#significance by fc & qval
	sig<-rep(0,length(q))
	
	sig[!is.na(fc) & fc>=fc_cutoff & !is.na(q) & q<q_cutoff]=1
	sig[!is.na(fc) & fc<=-fc_cutoff & !is.na(q) & q<q_cutoff]=-1
	sig[is.na(fc) | is.na(q) ]=NA

	mat.result<-cbind(fc,stat,p,q,sig)
	
	colnames(mat.result)<-c(paste0("Log2FC ",treat," vs ",ref),"NOISeq Stat:M,D,ranking","Q (1-prob)","Q (1-prob)",paste("Significance: Abs(Log2FC)>=",round(fc_cutoff,3)," ","Q<",q_cutoff,sep=""))
	
	result<-list()
	result$result=mat.result
	result$mydata=mydata
	result$res=res

	return(result)
}

######
#LM
######
lm_test <- function(mat,anno,design,fc_cutoff=1,q_cutoff=0.05,pmethod="Wald",qmethod="BH",core=core,treat,ref,independentfiltering=T,cookscutoff=T,na.rm=0){

	#Get the last variable
	design<-gsub(" ","",design)
	design.vars<-gsub("~","",design)
	design.vars<-unlist(strsplit(design.vars,"\\+"))
	
	#change treat and ref names
	treat.m=make.names(treat)
	ref.m=make.names(ref)

	#anno relevel column for last variable #only works for 1vs1 for the selected col
	#convert to factor
	anno[,design.vars[length(design.vars)]]<-factor(make.names(anno[,design.vars[length(design.vars)]]),levels=c(ref.m,treat.m))
	
	#new design
	design.mm<-data.frame(anno[,design.vars])
	colnames(design.mm)<-design.vars
	
	#remove NA and Inf
	is.na(mat) <- sapply(mat, is.infinite)
	mat[is.na(mat)]=na.rm	
	
	#lm test here
	#can be replaced by parallel implementation
	res<-t(apply(mat,1,function(x) {
		d<-data.frame(cbind(x,design.mm))
		lm.form<-as.formula(paste("x",design))
		result<-anova(lm(lm.form,d))
		c(result$`F value`[length(design.vars)],result$`Pr(>F)`[length(design.vars)])
	}))
	
	
	#convert to model matrix
	#design.mm<-model.matrix(as.formula(design),anno)	
	
	#remove empty col
	#design.mm<-design.mm[,apply(design.mm,2,sum)>0] 
	#print(design.mm)
	
	#remove NA
	#mat[is.na(mat)]=na.rm
	
	#lm test here
	#res<-t(apply(mat,1,function(x) {result<-anova(lm(x~design.mm));c(result$`F value`[1],result$`Pr(>F)`[1])}))
	
	#p,q,F
	p<-res[,2]
	q<-p.adjust(p,method=qmethod)
	stat<-res[,1]
	
	#calculate FC
	fc<-apply(mat,1,function(x){
		means<-as.vector(unlist(lapply(split(x,anno[,design.vars[length(design.vars)]]),mean)))
		if(means[1]>0) {
			if(means[2]==0) {
				0
			} else {
				log2(means[2]/means[1])
			}
		} else {
			log2((means[2]+0.001)/(means[1]+0.001))
		}
	})
	
	
	#significance by fc & qval
	sig<-rep(0,length(q))
	
	sig[!is.na(fc) & fc>=fc_cutoff & !is.na(q) & q<q_cutoff]=1
	sig[!is.na(fc) & fc<=-fc_cutoff & !is.na(q) & q<q_cutoff]=-1
	sig[is.na(fc) | is.na(q) ]=NA

	mat.result<-cbind(fc,stat,p,q,sig)
	
	colnames(mat.result)<-c(paste0("Log2FC ",treat," vs ",ref),"Stat:F","P",paste(qmethod,"P"),paste("Significance: ",treat," vs ",ref," Abs(Log2FC)>=",round(fc_cutoff,3)," ",qmethod, "P<",q_cutoff,sep=""))
	
	result<-list()
	result$result=mat.result

	return(result)
}


######
#Plots
######

volcano_plot_ggplot<-function(fc,q,sig,xlim=c(-7,7),ylim=c(0,30),xlab="Log2FC",ylab="-log10 P",main="Volcano Plot",fc_cutoff=args$fccutoff,q_cutoff=args$q_cutof) {
  
  q[is.na(q)]<-1 #remove NA
  
  fc<-as.numeric(unlist(fc))
  q<-as.numeric(unlist(q))
  
  #define color  
  cols <- c("Up" = "red", "Down" = "green","N.S."="grey")
  shs <- c("21" = 21, "24" = 24)
  #define col and shape
  
  shapes=rep(21,length(fc))
  shapes[abs(fc)>xlim[2]]=24
  shapes[-log10(q)>ylim[2]]=24
  
  sig.new<-rep("N.S.",length(sig))
  sig.new[sig==1]="Up"
  sig.new[sig==-1]="Down"
  
  
  #redefine cols and shs
  cols.sel<-cols[levels(as.factor(sig.new))]
  shs.sel<-shs[levels(as.factor(shapes))]
  
  #print(cols.sel)
  #print(shs.sel)
  
  #transform data
  fc[fc>xlim[2]]=xlim[2]
  fc[fc<xlim[1]]=xlim[1]
  
  q[-log10(q)>ylim[2]]=10^-ylim[2]
  
  #defined cols and shapes
  
  data<-data.frame( lfc=fc,q=-log10(q),sig=sig.new,shape=shapes)
  
  
  #plot
  vol <- ggplot(data, aes(x = lfc, y =q, fill = sig ,shape=factor(shape)))
  
  vol + ggtitle(label = main) +
    geom_point(size = 2, alpha = 1, na.rm = T, colour = "black") +
    scale_fill_manual(name="Color",values = cols) +
    scale_shape_manual(name="Shape",values = shs) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold",hjust = 0.5),axis.title=element_text(size=8)) + # change the legend
    guides(fill=guide_legend(override.aes = list(size = 3,  colour=cols.sel),title="Significance"),shape="none")+
    xlab(xlab) + # Change X-Axis label
    ylab(ylab) + # Change Y-Axis label
    ylim(ylim) +
    scale_x_continuous(breaks =seq(xlim[1],xlim[2],1),lim=xlim) +
    geom_hline(yintercept = -log10(as.numeric(q_cutoff)), colour="#990000", linetype="dashed") + #p cutoff
    geom_vline(xintercept = as.numeric(fc_cutoff), colour="#990000", linetype="dashed") + geom_vline(xintercept = -as.numeric(fc_cutoff), colour="#990000", linetype="dashed") + # fc cutoff line
    annotate("text",x=xlim[2]-0.5, y=ylim[2]-0.5, label=length(which(sig==1)),colour = "red",size = 5) + 
    annotate("text",x=xlim[1]+0.5, y=ylim[2]-0.5, label=length(which(sig==-1)),colour = "green",size = 5) + 
    annotate("text",x=xlim[2]-1, y=-log10(as.numeric(q_cutoff))+0.5, label=substring(main,28),colour = "red",size = 3) 
  
}

enhanced_volcano_plot <- function(gene, fc, q, sig, labels = NULL, 
                                  upcol = "red2", downcol = "blue2",
                                  xlab = "Log2FC", ylab = "-log10 P", main = "Volcano Plot", 
                                  fc_cutoff = args$fccutoff, q_cutoff = args$q_cutof,xlim=c(-10,10),ylim=c(0,30)){
    
  
	# Create a data frame using gene names, fold changes, and q-values
	fc <- as.numeric(unlist(fc))
	q <- as.numeric(unlist(q))
	q[is.na(q)] <- 1
	sig[is.na(sig)]=0

	gene <- make.names(as.character(unlist(gene)),unique = T) #changed here to allow duplicated gene names

	#remember the numbers
	up.number<-sum(sig==1)
	down.number<-sum(sig==-1)
	total.number<-length(sig)

	#decide xlim/ylim
	# Determine balanced x and y-axis limits (no features omitted)
	# y-axis upper limit
	if(min(q)>0) {
		max_pval <- -log10(min(q)) + 0.1
	} else {
		max_pval <- -log10(min(q[q>0])/10) + 0.1 #changed here for pval=0
	}
	
 	if(length(ylim)>1) {
		if(max_pval<max(ylim)) {
			ylim=c(0, max_pval)
		}
	} else {
		ylim=c(0, max_pval)
	} 
  
  
	# To get a symmetrical x-axis
	max_fc <- max(fc, na.rm = T)
	min_fc <- min(fc, na.rm = T)
	if (max_fc > abs(min_fc)){
		min_fc <- -max_fc
	}else{
		max_fc <- -min_fc
	}

	if(length(xlim)>1) {
		if(max_fc <max(xlim)) {
			xlim=c(min_fc, max_fc)
		}
	} else {
		xlim=c(min_fc, max_fc)
	}

	
	#filtered data
	newfc<-fc[abs(fc)<=max(xlim) & q>=10^-max(ylim)]
	newq<-q[abs(fc)<=max(xlim) & q>=10^-max(ylim)]
	newsig<-sig[abs(fc)<=max(xlim) & q>=10^-max(ylim)]
	newgene<-gene[abs(fc)<=max(xlim) & q>=10^-max(ylim)]
	
	
	#data used by volcano 
	df <- data.frame(gene_name = as.character(newgene), lfc = newfc, q = newq, stringsAsFactors = F)
	rownames(df) <- newgene
   

	# Create a named vector of custom colors

	col_scheme <- c("Up"=upcol, "Down"=downcol, "N.S."="grey")
	cols <- rep(col_scheme['N.S.'], length(newsig))
	cols[newsig == 1] <- col_scheme['Up']
	cols[newsig == -1] <- col_scheme['Down']
	names(cols)[cols == col_scheme['Up']] <- "Up"
	names(cols)[cols == col_scheme['Down']] <- "Down"
	names(cols)[cols == col_scheme['N.S.']] <- "N.S."
  
	# Label top 10 Up and down genes (based on FC) if user does not 
	# provide a vector of genes to label on volcano plot
	if (is.null(labels) == TRUE){
		tmp.df <- data.frame(newgene, newfc, cols)
		rownames(tmp.df) <- rownames(df)
		tmp.df <- tmp.df[tmp.df$cols != col_scheme['N.S.'],]
		tmp.df <- tmp.df[order(tmp.df$newfc),]
		labels <- c(rownames(tmp.df)[1:10], tail(rownames(tmp.df), n = 10))
	}
  

	# Render the volcano plot
	plt <- EnhancedVolcano(df, x = 'lfc', y = 'q', lab = df$gene_name,
						 pCutoff = q_cutoff, FCcutoff = fc_cutoff, 
						 gridlines.major = FALSE, gridlines.minor = FALSE, 
						 drawConnectors = F, legendLabSize = 12,
						 cutoffLineCol = "red", colAlpha = 0.75, 
						 cutoffLineType = "dashed", border = "full",
						 colCustom = cols, legendPosition = "right",
						 pointSize = 2, cutoffLineWidth = 0.4,
						 labFace = "plain", labSize=4,subtitle = main,
						 #ylim = c(0, max_pval), xlim = c(min_fc, max_fc),
						 ylim=ylim,xlim=xlim,
						 axisLabSize = 12, captionLabSize = 12, 
						 xlab = xlab, ylab = ylab, title = "", 
						 caption = paste0('Total = ', total.number, ' features'),
						 typeConnectors = "closed", legendIconSize = 2,
						 selectLab = labels, borderWidth = 1.5
						 )
  #changed drawConnectors = T,labSize=3
  
  # Add numbers of up and down DE genes to the plot
  
	plt <- plt + geom_text(x=min(xlim), y=max(ylim), label=down.number, 
						 col = col_scheme['Down'], size = 5)
	plt <- plt + geom_text(x=max(xlim), y=max(ylim), label=up.number, 
						 col = col_scheme['Up'], size = 5)
  
	print(plt) #generate plot
}

#to be implemented
ma_plot_ggplot<-function(m,a,sig,xlim=c(1,30),ylim=c(-7,7),xlab="A:Log2 Mean of Normalized Counts",ylab="M:Log2 Fold Change",main="MA Plot",m_cutoff=args$fc_cutof) {
  
  #m is lfc
  #a is mean
  
  #a[is.na(a)]<-1 #remove NA
  
  m<-as.numeric(unlist(m))
  a<-as.numeric(unlist(a))
  
  #define color  
  cols <- c("Up" = "red", "Down" = "blue","N.S."="grey")
  shs <- c("21" = 21, "24" = 24)
  #define col and shape
  
  shapes=rep(21,length(m))
  shapes[abs(m)>ylim[2]]=24
  shapes[a>xlim[2]]=24
  shapes[a<xlim[1]]=24
  
  sig.new<-rep("N.S.",length(sig))
  sig.new[sig==1]="Up"
  sig.new[sig==-1]="Down"
  
  
  #redefine cols and shs
  cols.sel<-cols[levels(as.factor(sig.new))]
  shs.sel<-shs[levels(as.factor(shapes))]
  
  #print(cols.sel)
  #print(shs.sel)
  
  #transform data
  a[a>xlim[2]]=xlim[2]
  a[a<xlim[1]]=xlim[1]
  
  m[m>ylim[2]]=ylim[2]
  m[m<ylim[1]]=ylim[1]
  
  #defined cols and shapes
  
  data<-data.frame( m=m,a=a,sig=sig.new,shape=shapes)
  
  
  #plot
  vol <- ggplot(data, aes(x = a, y =m, fill = sig ,shape=factor(shape)))
  
  vol <-vol + ggtitle(label = main) +
    geom_point(size = 2, alpha = 1, na.rm = T, colour = "black") +
    scale_fill_manual(name="Color",values = cols) +
    scale_shape_manual(name="Shape",values = shs) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold",hjust = 0.5),axis.title=element_text(size=8)) + # change the legend
    guides(fill=guide_legend(override.aes = list(size = 3,  colour=cols.sel),title="Significance"),shape="none")+
    xlab(xlab) + # Change X-Axis label
    ylab(ylab) + # Change Y-Axis label
    ylim(ylim) +
    scale_x_continuous(breaks =seq(xlim[1],xlim[2],1),lim=xlim) +
    #geom_hline(yintercept = -log10(as.numeric(a_cutoff)), colour="#990000", linetype="dashed") + #p cutoff
    geom_hline(yintercept = as.numeric(m_cutoff), colour="#990000", linetype="dashed") + geom_hline(yintercept = -as.numeric(m_cutoff), colour="#990000", linetype="dashed") + # a cutoff line
    annotate("text",x=xlim[2]-0.5, y=ylim[2]-0.5, label=length(which(sig==1)),colour = "red",size = 5) + 
    annotate("text",x=xlim[2]-0.5, y=ylim[1]+0.5, label=length(which(sig==-1)),colour = "blue",size = 5) 
	
	print(vol) #generate plot 
}

write_table_proper<-function(file,data,name="Gene") {
	data.df<-data.frame(name=rownames(data),data)
	names(data.df)<-c(name,colnames(data))

	write.table(data.df,file=file, row.names=FALSE,sep="\t",quote=F)

	#write.table(data.frame(name=rownames(data),data),file=file, row.names=FALSE,sep="\t",quote=F)
}

filter_de<-function(data,fc,q) {
	#filter the de result file for better volcano plot layout


}

#generate plots

#X 1. MA, 
#2. volcano
#3. hist for fc?
#4. ...


#####
#Statistical test
#####

#anno and config are changed in v0.6
data<-read.table(args$"in",header=T,row.names=1,sep="\t",check.names=F,flush=T,comment.char="",quote="")
config<-read.table(args$config,header=T,row.names=1,sep="\t",check.names=F,flush=T,colClasses="factor",comment.char="",quote="")

if(!is.na(args$anno)) {
	anno<-read.table(args$anno,header=T,row.names=1,sep="\t",check.names=F,comment.char="",quote="",flush=T)
	anno<-as.vector(unlist(anno[,1]))
} else {
	anno<-rownames(data)
}

#can be customized later
data.sel.rows<-filter_data(data,cutoff=args$filter)
data.sel<-data[data.sel.rows,]

anno.sel<-anno[data.sel.rows]


#save image
rdatafile=sub(".txt$",".rdata",args$out,perl=T)


if(args$pmethod == "DESeq2-Wald" | args$pmethod == "Wald") {
	data.sel.result<-deseq2_test(mat=data.sel,anno=config,design=args$formula,fc_cutoff=args$fccutoff,q_cutoff=args$qcutoff,pmethod=args$pmethod,qmethod=args$qmethod,treat=args$treat,ref=args$ref,independentfiltering=args$independentfiltering,cookscutoff=args$cookscutoff)
} else if (args$pmethod == "Lm") {
	data.sel.result<-lm_test(mat=data.sel,anno=config,design=args$formula,fc_cutoff=args$fccutoff,q_cutoff=args$qcutoff,pmethod=args$pmethod,qmethod=args$qmethod,treat=args$treat,ref=args$ref)
} else if (args$pmethod == "NOISeq" | args$pmethod == "NOIseq") {
	data.sel.result<-noiseq_test_norep(mat=data.sel,anno=config,design=args$formula,fc_cutoff=args$fccutoff,q_cutoff=args$qcutoff,treat=args$treat,ref=args$ref)
}

#write.table(data.sel.result$result,file=args$out,sep="\t",quote=F,col.names = NA)

write_table_proper(data.sel.result$result,file=args$out,"Feature")

save.image(file=rdatafile)


#####
#Plots
#####

if(args$plot) {

	
	#volcano plot
	vp_outfile_pdf=sub("\\.\\w+$","_volcanoplot.pdf",args$out,perl=T)

	#older version
	#pdf(vp_outfile_pdf,width=7, height=6)
	#volcano_plot_ggplot(fc=data.sel.result$result[,1],q=data.sel.result$result[,4],sig=data.sel.result$result[,5],xlab=colnames(data.sel.result$result)[1],ylab=paste("-Log10 ",colnames(data.sel.result$result)[4],sep=""),main=paste("Volcano Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),q_cutoff=args$qcutoff,fc_cutoff = args$fccutoff)
	#dev.off()


	#Volcano plot, pdf
	vp_outfile2_pdf=sub("\\.\\w+$","_volcanoplot2.pdf",args$out,perl=T)
	
	
	
	pdf(vp_outfile2_pdf,width=7, height=6)

	enhanced_volcano_plot(gene=anno.sel,fc=data.sel.result$result[,1],q=data.sel.result$result[,4],sig=data.sel.result$result[,5],xlab=colnames(data.sel.result$result)[1],ylab=paste("-Log10 ",colnames(data.sel.result$result)[4],sep=""),main=paste("Volcano Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),q_cutoff=args$qcutoff,fc_cutoff = args$fccutoff,upcol="red2",downcol="green2")

	dev.off()


	#Volcano plot, png
	vp_outfile_png=sub("\\.\\w+$","_volcanoplot.png",args$out,perl=T)

	CairoPNG(filename = vp_outfile_png,res = 300,width=2500, height=2200)
	
	data.sel.result.filtered<-filter_de(data.sel.result$result,fc=10,q=30)
	
	enhanced_volcano_plot(gene=anno.sel,fc=data.sel.result$result[,1],q=data.sel.result$result[,4],sig=data.sel.result$result[,5],xlab=colnames(data.sel.result$result)[1],ylab=paste("-Log10 ",colnames(data.sel.result$result)[4],sep=""),main=paste("Volcano Plot ","Significance: Abs(Log2FC)>=",round(args$fccutoff,2)," ",args$qmethod, "P<",args$qcutoff,sep=""),q_cutoff=args$qcutoff,fc_cutoff = args$fccutoff,xlim=c(-10,10),ylim=c(0,30))

	dev.off()

	vp_outfile_png2=sub("\\.\\w+$","_volcanoplot_nolim.png",args$out,perl=T)

	CairoPNG(filename = vp_outfile_png2,res = 300,width=2500, height=2200)

	enhanced_volcano_plot(gene=anno.sel,fc=data.sel.result$result[,1],q=data.sel.result$result[,4],sig=data.sel.result$result[,5],xlab=colnames(data.sel.result$result)[1],ylab=paste("-Log10 ",colnames(data.sel.result$result)[4],sep=""),main=paste("Volcano Plot ","Significance: Abs(Log2FC)>=",round(args$fccutoff,2)," ",args$qmethod, "P<",args$qcutoff,sep=""),q_cutoff=args$qcutoff,fc_cutoff = args$fccutoff,xlim="",ylim="")

	dev.off()

	#####
	#Plots
	#####

	#MA plot, pdf
	mp_outfile_pdf=sub("\\.\\w+$","_maplot.pdf",args$out,perl=T)

	pdf(mp_outfile_pdf,width=7, height=6)

	ma_plot_ggplot(m=data.sel.result$result[,1],a=log2(data.sel.result$res[,1]),sig=data.sel.result$result[,5],ylab=paste("M:",colnames(data.sel.result$result)[1],sep=""),main=paste("MA Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),m_cutoff = args$fccutoff)

	dev.off()

	#MA plot, jpg

	mp_outfile_png=sub("\\.\\w+$","_maplot.png",args$out,perl=T)

	CairoPNG(filename = mp_outfile_png,res = 300,width=2500, height=2200)

	ma_plot_ggplot(m=data.sel.result$result[,1],a=log2(data.sel.result$res[,1]),sig=data.sel.result$result[,5],ylab=paste("M:",colnames(data.sel.result$result)[1],sep=""),main=paste("MA Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),m_cutoff = args$fccutoff)

	dev.off()
	
	
	#mp_outfile_jpg=sub("\\.\\w+$","_maplot.jpg",args$out,perl=T)

	#CairoJPEG(filename = mp_outfile_jpg,res = 300)

	#ma_plot_ggplot(m=data.sel.result$result[,1],a=log2(data.sel.result$res[,1]),sig=data.sel.result$result[,5],ylab=paste("M:",colnames(data.sel.result$result)[1],sep=""),main=paste("MA Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),m_cutoff = args$fccutoff)
	#dev.off()
	
	#jpeg(mp_outfile_jpg,width=7, height=6, units="in", res=300)
	#ma_plot_ggplot(m=data.sel.result$result[,1],a=log2(data.sel.result$res[,1]),sig=data.sel.result$result[,5],ylab=paste("M:",colnames(data.sel.result$result)[1],sep=""),main=paste("MA Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),m_cutoff = args$fccutoff)
	#dev.off()
}
