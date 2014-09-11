install.packages("mvtnorm")
library("mvtnorm")
nts=9
nb.snps=10250
path=("~/Desktop/new_tutorial/")
setwd(path)
sum.stats=sapply(seq(1:nts),function(x){paste0(path,"out_eqtlbma_sumstats_tissue",x,".txt.gz",sep="")})


######Initialize Grid vals###
prior.grid=read.table("grid_weights.txt",header=F)
grid.vals=read.table("grid_phi2_oma2_with-configs.txt",header=F)
colnames(grid.vals)=c("phi2","omega2")
ugpo=list(NULL)
for(i in 1:length(grid.vals[,1])){##create a prior variance matrix for each set of grid values
  ugpo[[i]]=matrix(rep(grid.vals[i,2],nts),ncol=nts,nrow=nts)+diag(rep(grid.vals[i,1],nts))
}
########################
# read in summary statistics
#######################

sum.stat=list(NULL)
beta.hat=list(NULL)
sigma.hat=list(NULL)
for(i in 1:nts){
	df=read.table(sum.stats[i],header=T)
  sum.stat[[i]]=df
  beta.hat[[i]]=df$betahat.geno
  sigma.hat[[i]]=df$sebetahat.geno
}
  
############################################
#Cycle through each SNP gene pair according to user specified component ###
#################################
betas=rep(-1,1,by=0.01)##Values of beta to test
Id=diag(rep(1,nts))
post.var=list(NULL)
component.mu=list(NULL)
post.weight=matrix(NA,nrow=nb.snps,ncol=nrow(prior.grid))

for(p in 1:nb.snps){
  b.obs=unlist(lapply(beta.hat,function(x){return(x[p])}))###return beta  hats for gene-SNP pair across tissues
  sigma.obs=as.matrix(diag(unlist(lapply(sigma.hat,function(x){return(x[p])})),ncol=nts))###return sigma hats for gene-SNP pair across tissues
  pdata.l=(prior.grid*unlist(lapply(ugpo,function(x)(dmvnorm(b.obs,mean=rep(0,nts),sigma=x+sigma.obs))))) ##compute grid-specific pdatas
  pdata.total=sum(pdata.l)###compute pdata total
  post.ugp1=(lapply(ugpo,function(x) (solve(Id+x%*%solve(sigma.obs))%*%x)))##compute the posterior variance for each grid weight combo
  v.inv.beta.hat=solve(sigma.obs)%*%b.obs###helpful for posterior mean computation
  mu.gp1=lapply(post.ugp1,function(x)(x%*%v.inv.beta.hat))###Compute a list of expected posterior means vectors for each grid-weight combo
  post.weight[p,]=unlist(prior.grid/pdata.total)
  component.mu[[p]]=mu.gp1[[c]]##return user specificed mean vector
  post.var[[p]]=post.ugp1[[c]]###return user specified covariance matrix
 
}

################################  
  

	