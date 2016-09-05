# Generate data for paper tables and figs

rm(list = ls()) # Clean environment

if(!require(rgumlib)) {
  if(!require(devtools))
    install.packages("devtools",dependencies=TRUE)
  devtools::install_github("ppernot/rgumlib")
  library("rgumlib")
}
if(!require(numDeriv))
  install.packages("numDeriv",dependencies=TRUE) #  Numerical derivation
if(!require(SPIn))
  install.packages("SPIn",dependencies=TRUE) #  Shortest Probability Intervals

dir.create(path='results', 
           showWarnings = FALSE)

set.seed(1234) # Initialise la graine du RNG

models = c('ExpUnif','ExpNorm',
           'SinUnif','SinNorm',
           'TanhUnif','TanhUnif2')[1:6] 

for (model in models) {
  
  source(paste0('models/',model,'.R'))
  
  dx=(uMax-uMin)/(npt-1)
  res=matrix(0,ncol=6,nrow=npt)
  for ( i in 1:npt ) {
    x1.u = uMin +(i-1)*dx
    S=rgumlib::gumS1(nrun=nMC,adapt=FALSE,silent=TRUE,
                     fExpr=fExpr,x.mu=x1,x.u=x1.u,x.pdf=x1.pdf)
    res[i,1] = x1.u
    res[i,2] = sd(S$Y)
    res[i,3] = mean(S$Y)
    res[i,4] = rgumlib::vgSA(fExpr=fExpr,x.mu=x1,x.u=x1.u,X=S$X, Y=S$Y)$vg
    res[i,5] = diff(quantile(S$Y,probs=c(0.025,0.975),type=7))
    res[i,6] = diff(SPIn(S$Y)$spin)
  }
  save(fExpr,x1,x1.pdf,res,file=paste0('results/',model,'.Rda'))
}