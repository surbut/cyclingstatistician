**Understanding the Type Model: Sarah Urbut**

Briefly reintroducing the problem, suppose we have **G** genes upon whose expression we wish to assess in **S** tissues across **N** individuals, and that we also posess genetic information regarding their SNP genotype on these same **N** individuals. In our previous approach, we chose to jointly analyse the possibility that a particular SNP actived in *cis* as an eQTL across tissues types by enumerating all the possible configurations upon which a SNP could act and then computing a Bayes Factor $BF_{BMA}$ that quantified the likelihood of a SNP acting as an eQTL by averaging the $BF_{\gamma}$ likelihood overall possible configurations in which the SNP might be active. As the number of tissues considered in our previous analysis (Flutre *et al*, PLOS Genetics) becomes sufficiently large, it becomes computationally intractable to enumerate all the $2^{S}$ possible configurations $\gamma$. However, consider that the SNPs might actually fall into a set of 'types' in which given its membership in a particular type **K**, the SNP's predilection towards a particular configuration is generated by some underlying S-dimensional vector specific to the type, such that the SNP's *a priori* probability of being active in a particular tissue is independent of its activity in another tissue given it's type membership. That is, given an eQTL is of type k, its configuration $\mathbf \gamma$ is drawn from


$p_{k}( \mathbf \gamma| \mathbf q_{k})$ = $\prod {q_{ks}^{\gamma_{s}}(1−q_{ks})^{1−\gamma_{s}}}$

We then need only determine the number of types $K$, the vector of proportions in each type $\pi$, and the *S*-component vector of probabilities of being active in each tissue given membership in a type, $q_{k}$. Thus we have reduced the number of parameters to estimate from $2^{S}$ down to $K(S+1)$, because each type posesses the aforementioned vector of probabilities $q_{k}$ and the $\pi$ which corresponds to the probability a SNP is a member of each type. We hope this tutorial will help understanding this model by considering the generative capacity of such a framework.

As an overview,

1) For simplicity, we assume that (1-$\pi_{0}$) genes contain at most one eQTN and that this eQTN lies within 50 base pairs of the TSS. The following refer to SNPs in this category.

2) We then assign each SNP to a 'type' according to this probability vector $\mathbf \pi$, according the multinomial distribution with probability = $\mathbf \pi$ and size = 1 (because there is 1 type ID to sort into K potential types).

3) Conditional on membership in a type, we generate the *1* x *S* vector $q_{k,}$ containing the type-specific probability of being active in a particular tissue. Here, ach entry $q_{k,s}$ is drawn from the $\beta (1,1)$ distribution. The matrix *Q* contains each of the *K* *S* component vectors in its rows.

4) Now, conditional on this activity probability $q_{k}$ and type membership *k*, we generate a configuration of activity $\mathbf \gamma_{v}$ for each SNP. 

Now let's program!


```
## Warning: package 'MASS' was built under R version 3.0.2
```

```
## Installing package into '/Users/sarahurbut/Library/R/3.0/library'
## (as 'lib' is unspecified)
```

```
## Error: trying to use CRAN without setting a mirror
```

```
## Warning: package 'gtools' was built under R version 3.0.2
```

```r
nb.snps=100000
snps=seq(1,nb.snps)
# anno.snps=sapply(snps,function(x){
#   paste("rs",x,sep='')})
nb.genes=2000
##record gene.pos in terms of SNPs - e.g., a gene falls every 5 SNPS###
#gene.pos=seq(1,nb.snps,nb.snps/nb.genes)
n.inds=100
```
Here, we will use K=3 types and S=10 tissue types, but you can modify the situation for your needs. Similarly, for didactic purposes, we will assume that the majority of genes indeed contain an eQTL, thus $\pi_{0}$ is small (0.30) but you can generate randomly from the $\beta$ distribution or specify an alternative value. 

```r
K=3
S=10
#pi0=rbeta(1,1,1) ##prob gene has no eQTL
pi0=0.3
```
We should also consider covariates which might affect each tissue separately. Here, we consider only sex.

```r
set.seed(771)
gender=rbinom(n.inds,1,prob=0.5)
```
Furthermore we generate the vector of Minor Allele Frequencies (MAF) according to the Beta(1,1) distribution. We simulate genotype for each SNP according to this vector from the Binomial distribution.


```r
set.seed(707)
MAF=rbeta(nb.snps,1,1)
###Simulate genotypes in an n x v snps matrix such that each row corresponds to an individuals genotype vector##
n.geno=sapply(MAF,function(x)
  {rbinom(n=n.inds,prob=x,size=2)})
```

And herein lies the novelty of the method: we will simulate the prior probability of membership in a 'type' from the Dirichlet distribution. According to theses type proportions, we will then assign every SNP membership in a type.

```r
set.seed(771)
(pis=rdirichlet(1,alpha=rep(1,K)))
```

```
##       [,1]   [,2]   [,3]
## [1,] 0.376 0.2602 0.3638
```
Thus, about 38%, 26% and 36% of our SNPs are members of 'type' 1,2 or 3 respectively, which will drive their tissue specific expression patterns. 


```r
set.seed(302)
###simulate type identity according to pi proportion##
type=sapply(rep(0,nb.snps),function(x){
    assign.vec=rmultinom(n=1,size=1,prob=pis)
    which(assign.vec==1)}
    )
sum(type==1)/length(type)
```

```
## [1] 0.3738
```

```r
sum(type==2)/length(type)
```

```
## [1] 0.2597
```

```r
sum(type==3)/length(type)
```

```
## [1] 0.3665
```
We can check that the type proportion indeed match up with our vector of proportions $\mathbf \pi$.
Now, we wimulate type-specific vector $\mathbf q_{k}$, where **Q** is a **K** x **S** matrix indicating the probability of a SNP being an eQTL at each of the S tissues, given membership in type k.


```r
set.seed(111)
qmat=matrix(NA,nrow=K,ncol=S)
(q.mat=t(apply(qmat,1,function(x){
    x=rbeta(S,1,1)##generate predilection towards config from uniform
    })))
```

```
##        [,1]   [,2]    [,3]   [,4]   [,5]   [,6]   [,7]     [,8]   [,9]
## [1,] 0.4070 0.6296 0.62234 0.9893 0.5678 0.4442 0.9329 0.843797 0.8286
## [2,] 0.5689 0.6578 0.03247 0.3468 0.2126 0.9414 0.5342 0.640255 0.8837
## [3,] 0.3588 0.3643 0.42476 0.5634 0.3720 0.2748 0.9672 0.003283 0.4242
##       [,10]
## [1,] 0.6893
## [2,] 0.3579
## [3,] 0.9034
```
For example, for a SNP in type 1, the probability of being active in tissue 1 given it is the eQTL for the gene is 0.40, but the probability of being active in tissue 4 is close to 1 (0.98). We can see that the type of the SNP broadly dictates correlation patterns of activity between tissues.

Now, similar to my discussion of effect size and variance priors in the earlier tutorial, we simulate the covariance matrix $\mathbf \sigma_{GV}$ of the multivariate-normally-distributed vector $\mathbf beta_{s}$ according to the tissue-specific and tissue-average variance $\phi^{2}$ and $\omega^{2}$. We will draw the vector $\mathbf beta_{s}$ given a SNP's activity according the multivariate normal(0,$\mathbf \sigma_{GV}$) with the same covariance matrix for each gene. Similarly, for each gene, we simulate the covariance matrix of the errors between subgroup-expression levels as the *N* x *S* identity matrix *I*, assuming that the errors between tissues are uncorrelated.


```r
####Simulate covariance matrix of effect sizes, same for each gene##
oma2.plus.phi2 <- 0.8^2 # prior variance of subgroup SNP effect
oma2.over.oma2.plus.phi2 <- 3/4 # homogeneity
oma2 <- oma2.plus.phi2 * oma2.over.oma2.plus.phi2 # prior var of avg SNP effect
phi2 <- oma2.plus.phi2 * (1 - oma2.over.oma2.plus.phi2) # prior var of subgroup SNP effect, given avg effect
Sigma.beta.g <- matrix(rep(oma2, S^2),S, S) +
  diag(rep(phi2, S), S, S)
##Simulate covariance matrix of the errors, where we assume no covariance between tissues or individuals (i.e., matrix of multivariate normal S dimensional vectors###
###Key is assuming tissues and individuals are independent###
cov.err.S <- diag(S)
```
We also need to specify our matrix of covariates, using an inverse $\gamma$ prior on the tissue-specific variance from which we will draw our vector of tissue-specific effect size, $\beta_{c}$.


```r
set.seed(717)
###Create covariate matrix, using inverse gamma prior on variance ###
Sigma.beta.c <- diag(x=1 / rgamma(S, 2, 1),
                     nrow=S, ncol=S)
X.c <- cbind(rep(1, n.inds),gender)#use gender as a covariate
###TissueSpecific covariate effect##
B.c <- matrix(mvrnorm(n=2, mu=rep(0, S),
                      Sigma=Sigma.beta.c),
              nrow=2, ncol=S)
```
Now, for every gene, we first ask whether it is an 'eGene' - that is, is it affected by an eQTL. This is determined according to $\pi_{0}$, the proportion of null genes. Then, conditional on its differential expression, we will choose 1 'true eQTN', assuming 1 eQTN per gene, from set of all SNPS in *cis*, here defined as within 50 base pairs of the TSS. Given this SNPs type membership, we then generate a tissue-wide vecotr $\gamma$ dictating the behavior of this gene according to genotype at its eQTN from the type-specific probability vector $\q_{k}$.


```r
set.seed(500)
####Initialize DataFrames###
mean.exp=matrix(NA,nrow=nb.genes,ncol=S)
var.exp=matrix(NA,nrow=nb.genes,ncol=S)
cov.eQTL=matrix(NA,nrow=nb.genes,ncol=S)
truth=NULL

configurations=matrix(NA,nrow=nb.genes,ncol=S)
gene.type=NULL
for(i in 1:nb.genes){
  if(runif(1)<pi0){
    X=X.c
    B=B.c
  truth[i]=0
  configurations[i,] = rep(0,S)
  eqtn=0
  gene.type[i]=0
  }
else{
  
  start=50*(i-1)+1
  end=start+49
  cis=snps[start:end]##Specify the cis region
  eqtn=sample(cis,1)##sample the SNP from the chosen region, given that it contains an eQTN
  truth[i]=eqtn
  class=type[eqtn]##membership of the ith SNP
  gene.type[i]=class
  ##generate active tissues according to SNPs type membership#
  probs=q.mat[class,]
  config=sapply(probs,function(x){
  rbinom(1,1,prob=x)})
    configurations[i,] = config
    X.g <- matrix(n.geno[,eqtn], ncol=1)
    B.g <- matrix(mvrnorm(n=1, mu=rep(0, S),Sigma=Sigma.beta.g),nrow=1, ncol=S)
    B.g[1, which(config == 0)] <- 0
    E <- mvrnorm(n=n.inds, mu=rep(0,S),Sigma=cov.err.S)##draw each individuals ERROR vector from E_S(0,I)
    X <- cbind(X.c, X.g)
    B <- rbind(B.c, B.g)
}
  Y <- X %*% B + E
#print(eqtn)

  #write.table(Y, file=paste0(dir.name,i,"explevels.txt"), quote=F, sep="\t",
           # row.names=T, col.names=T)

mean.exp[i,]=apply(Y,2,mean)
var.exp[i,]=apply(Y,2,var)
if(eqtn==0){cov.eQTL[i,]=rep(0,S)}
else{cov.eQTL[i,]=apply(Y,2,function(x){cov(x,n.geno[,eqtn])})}
}
```
Note that we have the capacity to create, for every gene, an *N* by *S* matrix of expression levels across tissues for all individuals. We can then take the mean expression level for each gene across individuals for each tissue type, and record this in an *S* dimensional vector. Accordingly, we will have $M$ of these.

For each tissue, we can then look at boxplots of mean gene expression by type. We can compare these to the vector of configuration probabilities across types for a given tissue, i.e., *q[1:K,S]* with the expectation that for genes in a type in which the type-specific tissue probabilities $q_{k,s}$ is higher, there will be a greater variance in gene expression among the disase because the probability that $\gamma_{s}$ is 1 (and thus the eQTN is active given the gene contains an eQTN) is high. I plot this probability of tissue-specific activity above the type membership identifier.

```r
##########################
###order mean expression by gene.type###
gene.type=as.factor(gene.type)

plot.type=function(cov.eQTL,gene.type,tissue.types){
  par(mfrow=c(1,2),mar=c(5,5,2,2))
  for(i in tissue.types[1]:tissue.types[2]){
  #boxplot(mean.exp[,i]~gene.type,col=1:4,pch=1,xlab="GeneType",ylab="MeanExp",main=paste("tissue",i))
  #boxplot(var.exp[,i]~gene.type,col=1:4,pch=1,xlab="GeneType",ylab="MeanExp",main=paste("tissue",i))
  boxplot(abs(cov.eQTL[,i])~gene.type,ylim=c(-0.05,max(abs(cov.eQTL[,i]))),col=1:4,pch=1,xlab="GeneType",ylab="Cov(Gene,Geno)",main=paste("tissue",i))

  type.probs=round(q.mat[,i],2)
  text(x=2,y=0,type.probs[1],pos=1)
  text(x=3,y=0,type.probs[2],pos=1)
  text(x=4,y=0,type.probs[3],pos=1)
}
}

plot.type(cov.eQTL,gene.type,tissue.types=c(1,2))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) 

```r
plot.type(cov.eQTL,gene.type,tissue.types=c(3,4))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) 

```r
plot.type(cov.eQTL,gene.type,tissue.types=c(5,6))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-123.png) 

```r
plot.type(cov.eQTL,gene.type,tissue.types=c(7,8))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-124.png) 

```r
plot.type(cov.eQTL,gene.type,tissue.types=c(9,10))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-125.png) 

```r
install.packages("calibrate")
```

```
## Installing package into '/Users/sarahurbut/Library/R/3.0/library'
## (as 'lib' is unspecified)
```

```
## Error: trying to use CRAN without setting a mirror
```

```r
library(calibrate)
```

```
## Warning: package 'calibrate' was built under R version 3.0.2
```

```r
par(mfrow=c(2,5))
for(i in 1:S){
X=seq(1:K)
Y=sapply(seq(1:K),function(x){mean(abs(cov.eQTL[gene.type==x,i]))})
plot(X, Y,pch=1,col=2:4,type="p",ylim=c(0,0.5),xlim=c(0,3.5),cex=0.5,main="mean(abs(cov(gene.exp,geno)))",xlab="type",ylab="mean(abs(cov(gexp,geno)))")
library("calibrate")
textxy(X, Y, labs=round(q.mat[,i],2),cex=0.5,offset=.1)
       }
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-126.png) 

Neatly, let's consider tissue 3 for example. We can see that the covvariance in gene expression with genotype for genes whose eQTN are members of type 1 have a $q_{1,3}$ of 0.62. We would expect the majority of eQTNs from this type to exhibit active influence on gene expression in this tissue type. Correspondingly, we see that the covariance between gene expression and genotype is much higher for genes containing eQTN of this type than of genes containing eQTNs of type 2, for example. In type 2, the probability of being active in this particular tissue type given membership in type 2 is 0.03, and thus, the absolute value of the covariance between gene epxression and genotype is low.

**Applying CorMOTIF**

In order to apply this to cormotif, we need to look, gene be gene, at the expression levels for each genotype group. 
Our matrix will be $M$ by $3$ x $S$, because we will need to have a mean level of gene expression for each genotype type in all $S$ tissues.




However, the corMotif package aims to fit a model for differential expression between two groups, while eQTLs are really differential epxression between three groups.

This code pertains to our particular example, where I will show the flexibility of the model in the sense that within a type, not all SNPs mimic the same 'on' and off configuration in terms of differential expression. We can visualize the 'type-specific signature' similar to corMotif. However, CorMotif aims to capture type-specific differential expression across experiments, while we aim to capture type-specific genotype-associated-differential expression across tissue types. Though similar, our framework requires three levels of comparison for each subgroup, while CorMotif requires only two. I have modified their plotting framework to feature the reflect the proportion of genes in each type who show genotype-specific expression (i.e., the eQTL is active and $\gamma_{k,s}$ = 1). Thus darker bars indicate that the majority of Gene-SNP pairs in a type reflect an active eQTL relationship in this particular tissue, while lighter bars indicat the the majority of Gene-SNP pairs in a type are 'off' in the particular tissue and thus for most, $\gamma_{k,s}$ = 0). This is captured by the estimate of $q_{k,s}$, which reflects the proportion of SNPs in a type who are active as eQTLs in a particular tissue.


```r
ones=sum(gene.type==1)
  twos=sum(gene.type==2)
  threes=sum(gene.type==3)

freq.mat=apply(configurations,2,function(x){
  ones=sum(gene.type==1)
  twos=sum(gene.type==2)
  threes=sum(gene.type==3)
 c(sum(x[gene.type==1])/ones,
sum(x[gene.type==2])/twos,
sum(x[gene.type==3])/threes)
})



plot.eQTL=function (freq.mat, pis,title = "eQTL") 
{
    
  
  layout(matrix(1:2, ncol = 2))
    u <- 1:S##Nuber of Tissues
    v <- 1:K##Number of motifs
    image(u, v, t(freq.mat), col = gray(seq(from = 1, 
        to = 0, by = -0.1)), xlab = "Tissues", yaxt = "n", ylab = "Types", 
        main =paste(title, "pattern", sep = " "))
    axis(2, at = 1:length(v))
    
    for (i in 1:(length(u) + 1)) {
        abline(v = (i - 0.5))
    }
    for (i in 1:(length(v) + 1)) {
        abline(h = (i - 0.5))
    }
  Ng = 10000
  
  genecount = matrix(c(ones,twos,threes),nrow=1)
    NK = nrow(q.mat)
    plot(0, 0.7, pch = ".", xlim = c(0, 1.2), ylim = c(0.75, 
        NK + 0.25), frame.plot = FALSE, axes = FALSE, xlab = "No. of Gene-SNP pairs in each Class", 
        ylab = "", main = paste(title, "frequency", sep = " "))
    segments(0, 0.7,pis[1], 0.7)
    rect(0, 1:NK - 0.3, pis, 1:NK + 
        0.3, col = "dark grey")
    mtext(1:NK, at = 1:NK, side = 2, cex = 0.8)
    text(pis + 0.5, 1:NK, labels = paste(genecount,"(",round(pis,2),")",sep=""))
}


plot.eQTL(freq.mat,pis,title="eQTL")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 

This might be considered analogous to the framework of *Structure (Pritchard and Stephens)*, where we can think about the SNPs as the individuals belonging to an (unknown) number of subpopulations, here types. Accordingly, within a population, each individual exhibits a vector of alles at each locus drawn accoridng to the population-specific allele frequencies, with the probabilitiy of a given vector simply the product of the population-specific frequencies. Here, within a type, the 'SNPs' exhibit independent assignment of effect on gene-expression at each tissue, such that the probability of a given vector of activity $\gamma$ is the product of the type-specific probability at teach subgroup (or tissue).

