args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
load("essentials_allBRCA.RData")

integrand_e <- function(x,k) {dpois(k,x)}
integrand_m <- function(x,mean) {dnorm(x=mean,mean=x,sd=0.14)}

cat(paste("predicting for ",i,"\n",sep=""))
ptm <- proc.time()[3]
# prepare models and directories
system(intern=TRUE,command=paste('tar xf ',i,'.tar',sep=""))
system(command=paste('mkdir ./',i,'/toPredict',sep=""))

############### generate "missing" Var data ################
tempVar <- matrix(rep(".",length(samples)*3),nrow=length(samples),ncol=3)
colnames(tempVar) <- c("NAME:\tEXPR","M.GB","M.P")
rownames(tempVar) <- samples
eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/toPredict/samples_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
############### identify constitutive CpGs ################
IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
promoterVars <- promoterVars_template[1:length(IDs_promoter)]
promoter_CpGs <- template_promoter_CpGs[1:length(IDs_promoter)]
IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
geneBodyVars <- geneBodyVars_template[1:length(IDs_body)]
geneBody_CpGs <- template_body_CpGs[1:length(IDs_body)]
ncol = length(IDs_promoter) + length(IDs_body) + 1
############################################################

breaksEXPRESSION <- t(read.table(nrow=1,skip=1,file=paste("./",i,".result",sep="")))[,1]
breaksBODY <- t(read.table(nrow=1,skip=2,file=paste("./",i,".result",sep="")))[,1]
breaksPROMOTER <- t(read.table(nrow=1,skip=2,file=paste("./",i,".result",sep="")))[,1]

res_expr <- length(breaksEXPRESSION)-1
res_pr <- length(breaksPROMOTER)-1
res_gb <- length(breaksBODY)-1

all_labels_pr <- as.character(seq(1,res_pr,1))
all_labels_gb <- as.character(seq(1,res_gb,1))
all_labels_expr <- as.character(seq(1,res_expr,1))
promoter_samples <- matrix(ncol=res_pr,nrow=length(samples))
body_samples <- matrix(ncol=res_gb,nrow=length(samples))
expr_samples <- matrix(ncol=res_expr,nrow=length(samples))

epsilon_pr <- 1/(length(samples)*length(IDs_promoter))/res_pr
epsilon_gb <- 1/(length(samples)*length(IDs_body))/res_gb
epsilon_e <- 1/length(samples)/res_expr

tempS_samples <- matrix(ncol=ncol,nrow=length(samples))
for (current_sample in 1:length(samples)) {
	# expression
	read_count <- trunc(counts_BRCA[workingList_BRCA[i],samples[current_sample]])
	lambdas <- breaksEXPRESSION * factors_ls[current_sample]
	frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
	for (freq in 1:res_expr) {
		frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count,stop.on.error=FALSE)[1]
	}
	frequencies_expr <- unlist(frequencies_expr)
	if (length(which(frequencies_expr==0))==res_expr) frequencies_expr[length(frequencies_expr)] <- 1
	frequencies_expr <- frequencies_expr + epsilon_e
	frequencies_expr <- frequencies_expr/sum(frequencies_expr)
	
	# gene body
	cpg_list_gb <- NULL
	for (cpg in 1:length(geneBodyVars)) {
		miu <- mmatrix_pc[IDs_body[cpg],samples[current_sample]]
		if (!is.na(miu)) {
			frequencies_gb <- rep(0,res_gb)
			for (freq in 1:res_gb) {
			frequencies_gb[freq] <- integrate(integrand_m,lower=breaksBODY[freq],upper=breaksBODY[freq+1],mean=miu)$value
			}
			frequencies_gb <- unlist(frequencies_gb) + epsilon_gb
			frequencies_gb <- frequencies_gb/sum(frequencies_gb)
			cpg_list_gb[[cpg]] <- frequencies_gb
		} else cpg_list_gb[[cpg]] <- rep(1/res_gb,res_gb)
	}
	
	# promoter
	cpg_list_pr <- NULL
	for (cpg in 1:length(promoterVars)) {
		miu <- mmatrix_pc[IDs_promoter[cpg],samples[current_sample]]
		if (!is.na(miu)) {
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
			frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon_pr
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		} else cpg_list_pr[[cpg]] <- rep(1/res_pr,res_pr)
	}
	
	tempS_formated <- matrix(ncol=ncol,nrow=1)
	tempS_formated[1,1] <- paste('[1,',res_expr,']((',paste(frequencies_expr,sep="",collapse=","),'))',sep="",collapse="")
	cur_ncol <- 1
	for (element in 1:length(cpg_list_pr)) {
		tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="",collapse="")
		cur_ncol <- cur_ncol + 1
	}
	for (element in 1:length(cpg_list_gb)) {
		tempS_formated[1,cur_ncol+1] <- paste('[1,',res_gb,']((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep="",collapse="")
		cur_ncol <- cur_ncol + 1
	}
	tempS_samples[current_sample,] <- tempS_formated
}

tempFac <- tempS_samples
colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
rownames(tempFac) <- samples
eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/toPredict/samples_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))

# query the full AN model
string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/AN_model/all/ -l -n - ./',i,'/toPredict/samples_VarData.tab ./',i,'/toPredict/samples_FacData.tab',sep=""))
samples_G1_likelihoods <- as.numeric(substring(string[-1],17))

# query the full T model
string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/T_model/all/ -l -n - ./',i,'/toPredict/samples_VarData.tab ./',i,'/toPredict/samples_FacData.tab',sep=""))
samples_G2_likelihoods <- as.numeric(substring(string[-1],17))

# T posterior probability calculation
out <- cbind(samples,exp(-samples_G2_likelihoods) / (exp(-samples_G1_likelihoods) + exp(-samples_G2_likelihoods)),samples_G2_likelihoods,samples_G1_likelihoods)
colnames(out) <- c("sample_ID","posterior_G2","G2_mloglik","G1_mloglik")
eval(parse(text=paste('write.table(as.data.frame(out), col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE, file="./',i,'.predicted")',sep="")))

cat(paste("done predicting for ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
