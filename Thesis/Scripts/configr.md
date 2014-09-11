Computing the Posterior Probability Distribution on the Effect Size $\beta_{gp}$
========================================================

Following directly from the documentation on GitHub (configmodel), by maximum likelihood in each tissue separately, we can easily obtain the estimates of the standardized genotype effect sizes, $\hat{\beta}_{gp}$, and their standard errors recorded on the diagonal of an $R \times R$ matrix noted $\hat{V}_{gp} = Var(\hat{\beta}_{gp})$. 

This will be input by the user according to the summary stats function in eqtlbmaa_bf, as shown below with a sample from the tutorial. Here, the number of tissues was 9.


```r
install.packages("mvtnorm")
```

```
## Error: trying to use CRAN without setting a mirror
```

```r
library("mvtnorm")
```

```
## Warning: package 'mvtnorm' was built under R version 3.1.1
```

```r
nts=9
nb.snps=10250
path=("~/Desktop/new_tutorial/")
setwd(path)
sum.stats=sapply(seq(1:nts),function(x){paste0(path,"out_eqtlbma_sumstats_tissue",x,".txt.gz",sep="")})
sum.stat=list(NULL)
beta.hat=list(NULL)
sigma.hat=list(NULL)
for(i in 1:nts){
  df=read.table(sum.stats[i],header=T)
  sum.stat[[i]]=df
  beta.hat[[i]]=df$betahat.geno
  sigma.hat[[i]]=df$sebetahat.geno
}
```

Using each pair of tissues, we can also fill the off-diagonal elements of $\hat{V}_{gp}$. However, if the tissues arose from different studies, we may not have these covariance estimates available, and so we will specify that $\hat{V}_{gp}$ is diagonal for simplicity. If we now view $\hat{\beta}_{gp}$ and $\hat{V}_{gp}$ as *observations* (i.e. known), we can *forget* about the original data $X_p,X_c$ and $Y_g$, and  write a new **likelihood** (using only the sufficient statistics):

$\hat{\beta}_{gp} | \beta_{gp} \sim Norm_R(\beta_{gp}, \hat{V}_{gp})$

Let us imagine first that the prior on $\beta_{gp}$ is not a mixture but a single Normal: $\beta_{gp} \sim Norm_R(\mathbf{0}, U_{gp0})$.
As this prior is conjuguate to the ``likelihood'' above, the posterior simply is:

$\beta_{gp} | \hat{\beta}_{gp} \sim Norm_R({\mu}_{gp1}, U_{gp1})$

where:


$\mu_{gp1} = U_{gp1} (\hat{V}_{gp}^{-1} \hat{\beta}_{gp})$;
$U_{gp1} = (U_{gp0}^{-1} + \hat{V}_{gp}^{-1})^{-1}$.


In practice, we use a mixture for prior:

$\mathbf\beta_{gp}$ | $\lambda$, $\mathbf U_{gp0j}$ $\sim \sum_{l=1}^L \lambda_l \; Norm_R({0}, U_{gp0jl})$

for which the hyper-parameters are either fixed ($\mathbf U_{gp0j}$) or estimated ($\lambda$) as described using the full hierarchical model and the EM algorithm.
Moreover, the posteriors we seek are ${\beta}_{gp}$ | $\hat{\beta}_{gp}$, $\hat{V}_{gp}$, $\hat{\theta}$ (full marginal) or ${\beta}_{gp}$ | $\hat{\beta}_{gp}$, $\hat{V}_{gp}$, $\hat{\Theta}_{-\pi_0}$, $v_{gp}=1$ (conditional on being the eQTN in at least one tissue, averaging over all configurations and whole grid) or $\beta_{gp}$ | $\hat{\beta}_{gp}$, $\hat{V}_{gp}$, $\hat{\lambda}$, $\gamma_{gpr}=1$ (conditional on being an active eQTN in tissue $r$, averaging over all configurations in which this tissue is active and whole grid) or ${\beta}_{gp}$ | $\hat{\beta}_{gp}$, $\hat{V}_{gp}$, $\hat{\lambda}$, $c_{gpj}=1$ (conditional on being the eQTN and in configuration $j$, averaging over whole grid)$.

Below, we read in the prior grid weights and the corresponding combinations of heterogeneity, $\phi$ and $\omega$.


```r
######Initialize Grid vals###
prior.grid=read.table("grid_weights.txt",header=F)
grid.vals=read.table("grid_phi2_oma2_with-configs.txt",header=F)
colnames(grid.vals)=c("phi2","omega2")
prior.func=function(grid.vals,nts){
    ugpo=list(NULL)
    for(i in 1:length(grid.vals[,1])){##create a prior variance matrix for each set of grid values
  ugpo[[i]]=matrix(rep(grid.vals[i,2],nts),ncol=nts,nrow=nts)+diag(rep(grid.vals[i,1],nts))}
  return(ugpo)
  }
ugpo=prior.func(grid.vals,9)
```

Thus, the object *ugpo* now contains the 10 prior variance matrices for 9 tissues incroporating the values of $\phi$ and $\omega$ specified by our list of grid weights. Now, for each SNP-Gene Pair, we must consider the posterior distribution of effect sizes among the prespecified variety of grid sizes. Thus for every component, we will have l posterior distributions on $\beta_{gp}$ ~ $\hat{b}_{gp} \sim Norm_R({\mathbf \mu}_{gp1}, U_{gp1})$ and a corresponding set of posterior weights where each weight is proportion to the prior weight $\hat{lambda}$ and the probability of the data under all possible grid combinations.

$\sum_{l=1}^L \hat{\lambda}_l  \frac { Norm_R(\beta_{gp}|\mathbf{\mu_{gp1jl}}, U_{gp1jl})} {\sum_{l=1}^L \hat{\lambda}_l Norm_R(\hat{\beta}_{gp}|0, U_{gp0l} + \hat{V}_{gp})}$

Below, I create a function that, given a user-specified SNP-Gene pair, returns a vector of posterior weights and a list of posterior mean vectors and variance from this mixture distribution for every component.


```r
mixture.dist=function(p,nts){
  Id=diag(rep(1,nts))
 
  b.obs=unlist(lapply(beta.hat,function(x){return(x[p])}))###return beta  hats for gene-SNP pair across tissues
  sigma.obs=as.matrix(diag(unlist(lapply(sigma.hat,function(x){return(x[p])})),ncol=nts))###return sigma hats for gene-SNP pair across tissues
  pdata.l=(prior.grid*unlist(lapply(ugpo,function(x)(dmvnorm(b.obs,mean=rep(0,nts),sigma=x+sigma.obs))))) ##compute grid-specific pdatas
  pdata.total=sum(pdata.l)###compute pdata total
  post.ugp1=(lapply(ugpo,function(x) (solve(Id+x%*%solve(sigma.obs))%*%x)))##compute the posterior variance for each grid weight combo
  v.inv.beta.hat=solve(sigma.obs)%*%b.obs###helpful for posterior mean computation
  mu.gp1=lapply(post.ugp1,function(x)(x%*%v.inv.beta.hat))###Compute a list of expected posterior means vectors for each grid-weight combo
  post.weight=as.vector(unlist(prior.grid/pdata.total))
  post.weight=as.matrix(post.weight,nrow=1)
  return(list(component.means=mu.gp1,component.vars=post.ugp1,posterior.weights=unlist(post.weight)))}
```

Let's consider the posterior weights on $\beta_{gp=1}$ for example. We can also display the dataframe of grid weights corresponding to the covariance matrix from which they arose;


```r
mix=mixture.dist(p=1,nts=9)
post.weight=mix$posterior.weights
component.means=mix$component.means
component.vars=mix$component.vars
barplot(post.weight[,1],ylab="Posterior Component-Weight",xlab="ComponentID",main="Posterior Component Weights",col="blue")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

Then, we can choose which distribution to examine by calling considering the distributions of weights and choosing $\beta$ to arise from one of these according to the multinomial.


```r
set.seed(100)
component=which(rmultinom(1,1,prob=post.weight)==1)
print(component)
```

```
## [1] 7
```

Which shows our beta to be from component 7. Thus the distribution we expect for our $\beta$ is coming from the following distribution with mean vector and variance matrix:

```r
r=rmvnorm(100000,mean=component.means[[7]],component.vars[[7]])
component.means[[component]]
```

```
##            [,1]
##  [1,] -0.101562
##  [2,]  0.084140
##  [3,]  0.087395
##  [4,] -0.083268
##  [5,] -0.099651
##  [6,] -0.103920
##  [7,]  0.008849
##  [8,]  0.109624
##  [9,]  0.018479
```

```r
component.vars[[component]]
```

```
##           [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
##  [1,] 0.070811 0.004878 0.004703 0.004931 0.004792 0.004916 0.004700
##  [2,] 0.004878 0.071392 0.004739 0.004969 0.004829 0.004954 0.004736
##  [3,] 0.004703 0.004739 0.068664 0.004791 0.004656 0.004776 0.004566
##  [4,] 0.004931 0.004969 0.004791 0.072219 0.004881 0.005008 0.004787
##  [5,] 0.004792 0.004829 0.004656 0.004881 0.070049 0.004867 0.004652
##  [6,] 0.004916 0.004954 0.004776 0.005008 0.004867 0.071990 0.004773
##  [7,] 0.004700 0.004736 0.004566 0.004787 0.004652 0.004773 0.068612
##  [8,] 0.004986 0.005024 0.004844 0.005079 0.004936 0.005063 0.004841
##  [9,] 0.004861 0.004898 0.004723 0.004951 0.004812 0.004937 0.004719
##           [,8]     [,9]
##  [1,] 0.004986 0.004861
##  [2,] 0.005024 0.004898
##  [3,] 0.004844 0.004723
##  [4,] 0.005079 0.004951
##  [5,] 0.004936 0.004812
##  [6,] 0.005063 0.004937
##  [7,] 0.004841 0.004719
##  [8,] 0.073083 0.005007
##  [9,] 0.005007 0.071126
```

```r
par( mfrow = c( 3, 3 ), oma = c( 0, 0, 2, 0 ),mar=c(5,2,2,2) )

for(x in 1:nts){hist(r[,x],main="",freq=F,breaks=100,xlab=paste0("Tissue",x))}

title( "B|D for Component 7", outer = TRUE )
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

The vector $\mathbf r$ contains a set of randomly generated $\beta_{gp}$ drawn from the 7th component of our posterior distribution, according to the posterior weights and the corrresponding mean vector and posterior variance. 




