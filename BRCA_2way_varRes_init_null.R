args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1])
end <- as.numeric(args[2])
load("./essentials_BRCA.RData")
res_pr <- 50
res_gb <- 50
res_expr <- 50
epsilon <- 1e-300

tensor_product <- function(matrix1,matrix2) {
	if (is.matrix(matrix1) && is.matrix(matrix2)) result <- matrix(ncol=ncol(matrix1),nrow=ncol(matrix2),data=rep(0,ncol(matrix1)*ncol(matrix2))) else result <- matrix(ncol=length(matrix1),nrow=length(matrix1),data=rep(0,length(matrix1)*length(matrix2)))
	
	if (is.matrix(matrix1) && is.matrix(matrix2)) for (i in 1:nrow(matrix1)) {
		result <- result + matrix(ncol=ncol(matrix1),nrow=ncol(matrix2),byrow=TRUE,data=apply(expand.grid(matrix1[i,],matrix2[i,]), 1, prod))
	} else result <- result + matrix(nrow=length(matrix1),ncol=length(matrix2),byrow=TRUE,data=apply(expand.grid(matrix1,matrix2), 1, prod))
	
	#if (is.matrix(matrix1) && is.matrix(matrix2)) result <- result/nrow(matrix1)
	
	#for (i in 1:ncol(result)) result[,i] <- result[,i]/sum(result[,i])
	for (i in 1:nrow(result)) result[i,] <- result[i,]/sum(result[i,])
	#}
	return(result)
}

integrand_e <- function(x,k) {dpois(k,x)}
integrand_m <- function(x,mean) {dnorm(x=mean,mean=x,sd=0.14)}

geo_mean <- function(data) {
	log_data <- log(data)
	gm <- exp(mean(log_data[is.finite(log_data)]))
	return(gm)
}

M2Beta <- function (M) {
    return((2^M)/(2^M + 1))
}

Beta2M <- function (B) {
    return(log2(B/(1 - B)))
}

for (i in beg:end){
	cat(paste("doing ",i,"\n",sep=""))
	ptm <- proc.time()[3]
	system(command=paste('mkdir',i,sep=" "))
	system(command=paste('mkdir ./',i,'/AN_model',sep=""))
	system(command=paste('mkdir ./',i,'/AN_model/all',sep=""))
	system(command=paste('mkdir ./',i,'/T_model',sep=""))
	system(command=paste('mkdir ./',i,'/T_model/all',sep=""))
	system(command=paste('mkdir ./',i,'/full_model',sep=""))
	system(command=paste('mkdir ./',i,'/null',sep=""))
	system(command=paste('mkdir ./',i,'/null/AN_model',sep=""))
	system(command=paste('mkdir ./',i,'/null/T_model',sep=""))
	
	# generate "missing" Var data
	tempVar <- matrix(rep(".",length(ANs)*3),nrow=length(ANs),ncol=3)
	colnames(tempVar) <- c("NAME:\tEXPR","M.GB","M.P")
	rownames(tempVar) <- ANs
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/AN_model/all/AN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/null/AN_model/AN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
	rownames(tempVar) <- Ts
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/T_model/all/T_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/null/T_model/T_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	tempVar <- rbind(tempVar,tempVar)
	rownames(tempVar) <- c(Ts,ANs)
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/full_model/full_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	###########################################################
	############### identify constitutive CpGs ################
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	promoterVars <- promoterVars_template[1:length(IDs_promoter)]
	promoter_CpGs <- template_promoter_CpGs[1:length(IDs_promoter)]
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	geneBodyVars <- geneBodyVars_template[1:length(IDs_body)]
	geneBody_CpGs <- template_body_CpGs[1:length(IDs_body)]
	ncol = length(IDs_promoter) + length(IDs_body) + 1
	###########################################################
	############### define the binning scheme  ################
	all_labels_pr <- as.character(seq(1,res_pr,1))
	all_labels_gb <- as.character(seq(1,res_gb,1))
	all_labels_expr <- as.character(seq(1,res_expr,1))
	promoter_t <- matrix(ncol=res_pr,nrow=length(Ts))
	promoter_an <- matrix(ncol=res_pr,nrow=length(ANs))
	body_t <- matrix(ncol=res_gb,nrow=length(Ts))
	body_an <- matrix(ncol=res_gb,nrow=length(ANs))
	expr_t <- matrix(ncol=res_expr,nrow=length(Ts))
	expr_an <- matrix(ncol=res_expr,nrow=length(ANs))
	pseudo_counts_pr <- matrix(ncol=res_pr,nrow=res_expr,data=rep(1,res_expr*res_pr))
	pseudo_counts_gb <- matrix(ncol=res_gb,nrow=res_expr,data=rep(1,res_expr*res_gb))
	
	# gene body
	x <- mmatrix[IDs_body,c(Ts,ANs)]
	
	density <- density(x,bw=0.14,from=-7,to=7,n=2801)
	density$y <- density$y/sum(density$y)
	density$y <- cumsum(density$y)
	breaks <- NULL
	noBreaks <- res_gb-1
	for (j in 1:noBreaks) { breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
	breaksBODY <- sort(c(-7.01,breaks,7.01))
	
	# promoter
	x <- mmatrix[IDs_promoter,c(Ts,ANs)]
	
	density <- density(x,bw=0.14,from=-7,to=7,n=2801)
	density$y <- density$y/sum(density$y)
	density$y <- cumsum(density$y)
	breaks <- NULL
	noBreaks <- res_pr-1
	for (j in 1:noBreaks) { breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
	breaksPROMOTER <- sort(c(-7.01,breaks,7.01))
	
	# expression
	temp <- counts_BRCA_plusOne[workingList_BRCA[i],c(Ts,ANs)]
	temp <- as.data.frame(temp)
	colnames(temp) <- "EXPRESSION"
	tempAN <- matrix(ncol=2)
	colnames(tempAN) <- c("cpm","density")
	for (j in 1:nrow(temp)) {
		lambda <- temp[j,1]
		X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
		current <- factors_ls[c(which_Ts,which_ANs)[j]]
		tempAN <- rbind(tempAN,cbind(X/current,dpois(X,lambda=lambda)*current))
	}
	tempAN <- as.data.frame(tempAN[-1,],)
	tempAN <- tempAN[order(tempAN$cpm),]
	tempAN[,3] <- cumsum(tempAN[,2])
	tempAN[,3] <- tempAN[,3]/max(tempAN[,3])
	breaks <- NULL
	noBreaks <- res_expr-1
	for (j in 1:noBreaks) { breaks <- c (breaks, tempAN[which(tempAN[,3] >= j*(1/(1+noBreaks))),1][1])}
	breaksEXPRESSION <- c(0,breaks,Inf)
	########################################################
	# dynamic generation of model specification files here #
	########################################################
	####################### state map ######################
	stateMaps <- file(paste("./",i,"/stateMaps.txt",sep=""),"w")
	exprMap <- paste("NAME:\texprMap\nSYMBOLS:\t",paste(seq(1,res_expr,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_expr,1),collapse=" "),"; *=",paste(seq(1,res_expr,1),collapse=" "),";\n\n",collapse="",sep="")
	prMap <- paste("NAME:\tprMap\nSYMBOLS:\t",paste(seq(1,res_pr,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_pr,1),collapse=" "),"; *=",paste(seq(1,res_pr,1),collapse=" "),";\n\n",collapse="",sep="")
	gbMap <- paste("NAME:\tgbMap\nSYMBOLS:\t",paste(seq(1,res_gb,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_gb,1),collapse=" "),"; *=",paste(seq(1,res_gb,1),collapse=" "),";\n",collapse="",sep="")
	cat(exprMap,prMap,gbMap,file=stateMaps)
	close(stateMaps)
	########################################################
	####################### variables ######################
	variables <- file(paste("./",i,"/variables.txt",sep=""),"w")
	cat("STATE_MAP_NAME:\texprMap\nVAR_NAMES:\tEXPR\n\n",sep="",file=variables)
	cat("STATE_MAP_NAME:\tprMap\nVAR_NAMES:\tM.P\n\n",sep="",file=variables)
	cat("STATE_MAP_NAME:\tgbMap\nVAR_NAMES:\tM.GB\n",sep="",file=variables)
	close(variables)
	########################################################
	##################### factor graph #####################
	factorGraph <- file(paste("./",i,"/factorGraph.txt",sep=""),"w")
	cat("NAME:\tEXPR.likelihood\nNB1:\tEXPR\nPOT:\tpot_EXPR.likelihood\n",file=factorGraph)
	cat("\nNAME:\tEXPR.prior\nNB1:\tEXPR\nPOT:\tpot_EXPR.prior\n",file=factorGraph)
	cat(paste("\nNAME:\tEXPR.M.GB\nNB1:\tEXPR\nNB2:\tM.GB\nPOT:\tpot_EXPR.M.GB\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\tEXPR.M.P\nNB1:\tEXPR\nNB2:\tM.P\nPOT:\tpot_EXPR.M.P\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\t",promoter_CpGs,"\nNB1:\tM.P\nPOT:\tpot_",promoter_CpGs,"\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\t",geneBody_CpGs,"\nNB1:\tM.GB\nPOT:\tpot_",geneBody_CpGs,"\n",sep="",collapse=""),file=factorGraph)
	close(factorGraph)
	
	system(command=paste('cp ./',i,'/*.txt ./',i,'/AN_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/T_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/full_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/AN_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/T_model',sep=""))
	###########################################################################
	############################## AN model ###################################
	## full AN model developed from here, to obtain likelihoods of ANs ########
	
	# generate FacData for full set of ANs
	tempS_AN <- matrix(ncol=ncol,nrow=length(ANs))
	for (current_sample in 1:length(ANs)) {
		# expression
		read_count <- trunc(counts_BRCA_plusOne[workingList_BRCA[i],ANs[current_sample]])
		lambdas <- breaksEXPRESSION * factors_ls[which_ANs[current_sample]]
		frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
		for (freq in 1:res_expr) {
			frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count)[1]
		}
		frequencies_expr <- unlist(frequencies_expr)
		if (length(which(frequencies_expr==0))==5) frequencies_expr[length(frequencies_expr)] <- 1
		frequencies_expr <- frequencies_expr + epsilon
		
		# gene body
		cpg_list_gb <- NULL
		for (cpg in 1:length(geneBodyVars)) {
			miu <- mmatrix[IDs_body[cpg],ANs[current_sample]]
			frequencies_gb <- rep(0,res_gb)
			for (freq in 1:res_gb) {
				frequencies_gb[freq] <- integrate(integrand_m,lower=breaksBODY[freq],upper=breaksBODY[freq+1],mean=miu)$value
			}
			frequencies_gb <- unlist(frequencies_gb) + epsilon
			cpg_list_gb[[cpg]] <- frequencies_gb
		}
		
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix[IDs_promoter[cpg],ANs[current_sample]]
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon
			cpg_list_pr[[cpg]] <- frequencies_pr
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
		tempS_AN[current_sample,] <- tempS_formated
		
		#start precomputing correct initialization of parameters
		promoter_an[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
		body_an[current_sample,] <- apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean))
		expr_an[current_sample,] <- frequencies_expr
	}
	
	# build and the full model with AN samples
	# precompute correct initialization of parameters for AN-only model
	prior_pr <- apply(promoter_an,2,mean)/sum(apply(promoter_an,2,mean))
	prior_gb <- apply(body_an,2,mean)/sum(apply(body_an,2,mean))
	prior_expr <- apply(expr_an,2,mean)/sum(apply(expr_an,2,mean))
	
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_gb,collapse=",")
	geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(body_an,expr_an)
	expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(promoter_an,expr_an)
	expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
	
	potentials <- file(paste("./",i,"/AN_model/all/factorPotentials.txt",sep=""),"w")
	cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_AN
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
	rownames(tempFac) <- ANs
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/AN_model/all/AN_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	#  query the full model with AN samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/AN_model/all/ -l -n - ./',i,'/AN_model/all/AN_VarData.tab ./',i,'/AN_model/all/AN_FacData.tab',sep=""))
	ANs_AN_likelihoods <- as.numeric(substring(string[-1],17))
	##########################################################################
	############################## T model ###################################
	### full T model developed from here, to obtain likelihoods of Ts #######

	# generate FacData for full set of Ts
	tempS_T <- matrix(ncol=ncol,nrow=length(Ts))
	for (current_sample in 1:length(Ts)) {
		# expression
		read_count <- trunc(counts_BRCA_plusOne[workingList_BRCA[i],Ts[current_sample]])
		lambdas <- breaksEXPRESSION * factors_ls[which_Ts[current_sample]]
		frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
		for (freq in 1:res_expr) {
			frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count)[1]
		}
		frequencies_expr <- unlist(frequencies_expr)
		if (length(which(frequencies_expr==0))==5) frequencies_expr[length(frequencies_expr)] <- 1
		frequencies_expr <- frequencies_expr + epsilon
		
		# gene body
		cpg_list_gb <- NULL
		for (cpg in 1:length(geneBodyVars)) {
			miu <- mmatrix[IDs_body[cpg],Ts[current_sample]]
			frequencies_gb <- rep(0,res_gb)
			for (freq in 1:res_gb) {
				frequencies_gb[freq] <- integrate(integrand_m,lower=breaksBODY[freq],upper=breaksBODY[freq+1],mean=miu)$value
			}
			frequencies_gb <- unlist(frequencies_gb) + epsilon
			cpg_list_gb[[cpg]] <- frequencies_gb
		}
		
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix[IDs_promoter[cpg],Ts[current_sample]]
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon
			cpg_list_pr[[cpg]] <- frequencies_pr
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
		tempS_T[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_t[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
		body_t[current_sample,] <- apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean))
		expr_t[current_sample,] <- frequencies_expr
	}
	
	# build and the full model with AN samples
	# precompute correct initialization of parameters for AN-only model
	prior_pr <- apply(promoter_t,2,mean)/sum(apply(promoter_t,2,mean))
	prior_gb <- apply(body_t,2,mean)/sum(apply(body_t,2,mean))
	prior_expr <- apply(expr_t,2,mean)/sum(apply(expr_t,2,mean))
	
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_gb,collapse=",")
	geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(body_t,expr_t)
	expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(promoter_t,expr_t)
	expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
	
	potentials <- file(paste("./",i,"/T_model/all/factorPotentials.txt",sep=""),"w")
	cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_T
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
	rownames(tempFac) <- Ts
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/T_model/all/T_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	#  query the full model with T samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/T_model/all/ -l -n - ./',i,'/T_model/all/T_VarData.tab ./',i,'/T_model/all/T_FacData.tab',sep=""))
	Ts_T_likelihoods <- as.numeric(substring(string[-1],17))
	##########################################################################
	
	###########################################################################
	######################## All data model ###################################
	## Full model developed from here, to obtain likelihoods of Ts and ANs ####
	# precompute correct initialization of parameters for AN-only model
	promoter_all <- rbind(promoter_t,promoter_an)
	body_all <- rbind(body_t,body_an)
	expr_all <- rbind(expr_t,expr_an)
	
	prior_pr <- apply(promoter_all,2,mean)
	prior_gb <- apply(body_all,2,mean)
	prior_expr <- apply(expr_all,2,mean)
	
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_gb,collapse=",")
	geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(body_all,expr_all)
	expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(promoter_all,expr_all)
	expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
	
	potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
	cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
	close(potentials)
	
	tempFac <- rbind(tempS_T,tempS_AN)
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
	rownames(tempFac) <- c(Ts,ANs)
	eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with T and AN samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/full_model/ -l -n - ./',i,'/full_model/full_VarData.tab ./',i,'/full_model/full_FacData.tab',sep=""))
	allData_full_likelihoods <- as.numeric(substring(string[-1],17))
	###########################################################################################
	
	###########################################################################
	######################## D calculation ###################################
	D <- 2*(sum(allData_full_likelihoods) - (sum(ANs_AN_likelihoods)+sum(Ts_T_likelihoods)))
	###########################################################################
	################# P val calculation using null distr. #####################
	nruns <- 100
	Ds <- vector(length=nruns,mode="numeric")
	for (run in 1:nruns) {
		cur <- sample(x=1:(length(Ts)+length(ANs)),size=length(Ts),replace=FALSE)
		# Ts
		tempFac_T <- tempFac[cur,]
		rownames(tempFac_T) <- Ts
		eval(parse(text = paste('write.table(', paste('tempFac_T,file = "./',i,'/null/T_model/T_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		prior_pr <- apply(promoter_all[cur,],2,mean)
		prior_gb <- apply(body_all[cur,],2,mean)
		prior_expr <- apply(expr_all[cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_all[cur,],expr_all[cur,])
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_all[cur,],expr_all[cur,])
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/null/T_model/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/T_model/ -l -n - ./',i,'/T_model/all/T_VarData.tab ./',i,'/null/T_model/T_FacData.tab',sep=""))
		Ts_T_likelihoods <- as.numeric(substring(string[-1],17))
		
		# ANs
		tempFac_AN <- tempFac[-cur,]
		rownames(tempFac_AN) <- ANs
		eval(parse(text = paste('write.table(', paste('tempFac_AN,file = "./',i,'/null/AN_model/AN_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		prior_pr <- apply(promoter_all[-cur,],2,mean)
		prior_gb <- apply(body_all[-cur,],2,mean)
		prior_expr <- apply(expr_all[-cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_all[-cur,],expr_all[-cur,])
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_all[-cur,],expr_all[-cur,])
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/null/AN_model/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/AN_model/ -l -n - ./',i,'/AN_model/all/AN_VarData.tab ./',i,'/null/AN_model/AN_FacData.tab',sep=""))
		ANs_AN_likelihoods <- as.numeric(substring(string[-1],17))
		
		Ds[run] <- 2*(sum(allData_full_likelihoods) - (sum(ANs_AN_likelihoods)+sum(Ts_T_likelihoods)))
	}
	pval_zscore <- 1-pnorm(D,mean=mean(Ds),sd=sd(Ds))
	zscore <- (D - mean(Ds)) / sd(Ds)
	###########################################################################################
	
	eval(parse(text=paste('write.table(x=t(c(pval_zscore,D,mean(Ds),sd(Ds),zscore)), col.names=FALSE, row.names=FALSE, append=FALSE, file="./',i,'.result")',sep="")))
	cat(paste("done ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
