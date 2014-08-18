args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1])
end <- as.numeric(args[2])
load("./essentials_BRCA.RData")
library(aws)
res_pr <- 50
epsilon <- 1e-300

integrand_m <- function(x,mean) {dnorm(x=mean,mean=x,sd=0.14)}

geo_mean <- function(data) {
	log_data <- log(data)
	gm <- exp(mean(log_data[is.finite(log_data)]))
	return(gm)
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
	
	IDs_promoter <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	promoterVars <- geneBodyVars_template[1:length(IDs_promoter)]
	promoter_CpGs <- template_body_CpGs[1:length(IDs_promoter)]
	# generate "missing" Var data
	ncol = length(IDs_promoter)
	tempVar <- matrix(rep(".",length(ANs)*1),nrow=length(ANs),ncol=1)
	colnames(tempVar) <- "NAME:\tM.P"
	rownames(tempVar) <- ANs
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/AN_model/all/AN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
	rownames(tempVar) <- Ts
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/T_model/all/T_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	tempVar <- rbind(tempVar,tempVar)
	rownames(tempVar) <- c(Ts,ANs)
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/full_model/full_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	###########################################################
	############### binning scheme defined here ###############
	# promoter
	x <- as.vector(mmatrix[IDs_promoter,c(Ts,ANs)])
	
	density <- density(x,bw=0.14,from=-7,to=7,n=2801)
	density$y <- density$y/sum(density$y)
	density$y <- cumsum(density$y)
	breaks <- NULL
	noBreaks <- res_pr-1
	for (j in 1:noBreaks) { breaks <- c(breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
	breaksPROMOTER <- sort(c(-7.01,breaks,7.01))
	
	all_labels_pr <- as.character(seq(1,res_pr,1))
	promoter_t <- matrix(ncol=res_pr,nrow=length(Ts))
	promoter_an <- matrix(ncol=res_pr,nrow=length(ANs))
	########################################################
	# dynamic generation of model specification files here #
	########################################################
	stateMaps <- file(paste("./",i,"/stateMaps.txt",sep=""),"w")
	prMap <- paste("NAME:\tprMap\nSYMBOLS:\t",paste(seq(1,res_pr,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_pr,1),collapse=" "),"; *=",paste(seq(1,res_pr,1),collapse=" "),";\n\n",collapse="",sep="")
	cat(prMap,file=stateMaps)
	close(stateMaps)
	########################################################
	####################### variables ######################
	variables <- file(paste("./",i,"/variables.txt",sep=""),"w")
	cat("STATE_MAP_NAME:\tprMap\nVAR_NAMES:\tM.P\n\n",sep="",file=variables)
	close(variables)
	########################################################
	##################### factor graph #####################
	factorGraph <- file(paste("./",i,"/factorGraph.txt",sep=""),"w")
	cat(paste("\nNAME:\tprior.M.P\nNB1:\tM.P\nPOT:\tpot_P.M.prior\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\t",promoter_CpGs,"\nNB1:\tM.P\nPOT:\tpot_",promoterVars,"\n",sep="",collapse=""),file=factorGraph)
	close(factorGraph)
	
	system(command=paste('cp ./',i,'/*.txt ./',i,'/AN_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/T_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/full_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/AN_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/T_model',sep=""))
	
	promoter_CpGs[1] <- paste("NAME:\t",promoter_CpGs[1],sep="")
	###########################################################################
	############################## AN model ###################################
	## full AN model developed from here, to obtain likelihoods of ANs ########
	
	# generate FacData for full set of ANs
	tempS_AN <- matrix(ncol=ncol,nrow=length(ANs))
	for (current_sample in 1:length(ANs)) {
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix[IDs_promoter[cpg],ANs[current_sample]]
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		cur_ncol <- 0
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_AN[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_an[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
	}
	# precompute correct initialization of parameters for AN-only model
	prior_pr <- kernsm(apply(promoter_an,2,mean),h=2)
	prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/AN_model/all/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_AN
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- ANs
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/AN_model/all/AN_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with AN samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/AN_model/all/ -l -n - ./',i,'/AN_model/all/AN_VarData.tab ./',i,'/AN_model/all/AN_FacData.tab',sep=""))
	ANs_AN_likelihoods <- as.numeric(substring(string[-1],17))
	##########################################################################
	############################## T model ###################################
	### full T model developed from here, to obtain likelihoods of Ts #######
	tempS_T <- matrix(ncol=ncol,nrow=length(Ts))
	for (current_sample in 1:length(Ts)) {
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix[IDs_promoter[cpg],Ts[current_sample]]
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		cur_ncol <- 0
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_T[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_t[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
	}
	# precompute correct initialization of parameters for T-only model
	prior_pr <- kernsm(apply(promoter_t,2,mean),h=2)
	prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/T_model/all/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_T
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- Ts
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/T_model/all/T_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with T samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/T_model/all/ -l -n - ./',i,'/T_model/all/T_VarData.tab ./',i,'/T_model/all/T_FacData.tab',sep=""))
	Ts_T_likelihoods <- as.numeric(substring(string[-1],17))
	###########################################################################
	######################## All data model ###################################
	## Full model developed from here, to obtain likelihoods of Ts and ANs ####
	# precompute correct initialization of parameters for full model
	promoter_all <- rbind(promoter_t,promoter_an)
	
	prior_pr <- kernsm(apply(promoter_all,2,mean),h=2)
	prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)

	tempFac <- rbind(tempS_T,tempS_AN)
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- c(Ts,ANs)
	eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# build and query the full model with T and AN samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/full_model/ -l -n - ./',i,'/full_model/full_VarData.tab ./',i,'/full_model/full_FacData.tab',sep=""))
	allData_full_likelihoods <- as.numeric(substring(string[-1],17))
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
		tempFac_T <- as.data.frame(tempFac[cur,])
		colnames(tempFac_T) <- promoter_CpGs[1:length(promoterVars)]
		rownames(tempFac_T) <- Ts
		eval(parse(text = paste('write.table(', paste('tempFac_T,file = "./',i,'/null/T_model/T_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		prior_pr <- kernsm(apply(promoter_all[cur,],2,mean),h=2)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
		potentials <- file(paste("./",i,"/null/T_model/factorPotentials.txt",sep=""),"w")
		cat(promoterPots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/T_model/ -l -n - ./',i,'/T_model/all/T_VarData.tab ./',i,'/null/T_model/T_FacData.tab',sep=""))
		Ts_T_likelihoods <- as.numeric(substring(string[-1],17))
		
		# ANs
		tempFac_AN <- as.data.frame(tempFac[-cur,])
		colnames(tempFac_AN) <- promoter_CpGs[1:length(promoterVars)]
		rownames(tempFac_AN) <- ANs
		eval(parse(text = paste('write.table(', paste('tempFac_AN,file = "./',i,'/null/AN_model/AN_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		prior_pr <- kernsm(apply(promoter_all[-cur,],2,mean),h=2)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
		potentials <- file(paste("./",i,"/null/AN_model/factorPotentials.txt",sep=""),"w")
		cat(promoterPots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/AN_model/ -l -n - ./',i,'/AN_model/all/AN_VarData.tab ./',i,'/null/AN_model/AN_FacData.tab',sep=""))
		ANs_AN_likelihoods <- as.numeric(substring(string[-1],17))
		# Ds calculation
		Ds[run] <- 2*(sum(allData_full_likelihoods) - (sum(ANs_AN_likelihoods)+sum(Ts_T_likelihoods)))
	}
	if (D != 0) pval_zscore <- 1-pnorm(D,mean=mean(Ds),sd=sd(Ds)) else pval_zscore <- 1
	if (sd(Ds) != 0) zscore <- (D - mean(Ds)) / sd(Ds) else zscore <- 0
	###########################################################################################
	eval(parse(text=paste('write.table(x=t(c(pval_zscore,D,mean(Ds),sd(Ds),zscore)), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	cat(paste("done ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
