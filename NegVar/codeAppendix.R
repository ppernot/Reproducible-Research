vgSA = function(yFunc,X,Y) {
  # Compute VG indices 
  # - yFunc : model function 
  # - X     : matrix of inputs values 
  # - Y     : vector of output values
  if(!require(numDeriv))
    install.packages("numDeriv",dependencies=TRUE) #  Numerical derivation
  
  # Compute derivatives 
  fD = apply(X, 1,                
             function(x) {                 
               fun = function(x) 
                 do.call(yFunc,as.list(x))
               numDeriv::grad(fun,x)               
             }   
  )
  tabDeriv = matrix(fD,ncol=ncol(X),
                    byrow=TRUE) 
  
  # Center samples 
  Xc = scale(X,center=TRUE,scale=FALSE)
  Yc = scale(Y,center=TRUE,scale=FALSE)
  
  # Build matrix of outputs conforming with Xc
  Yc = matrix(Yc,ncol=ncol(Xc),nrow=nrow(Xc), 
              byrow=FALSE)
  
  # Compute VGs
  vg = colMeans(Yc*tabDeriv*Xc)/var(Y)
  
  return(vg) 
}

# Monte Carlo sample size 
nMC = 10^6

# Model function 
yFunc = function(x1) exp(-abs(x1))

# Define uncertain variable(s) 
ux = 2/3^0.5  # Standard uncertainty of X
a = ux*3^0.5  # Bounds of uniform distrib. 

# Sample uncertain variable(s) 
X = matrix(runif(nMC, min=-a, max=a),ncol=1) 

# Ensure the match between sample
# columns and yFunc argument(s)
colnames(X) = 'x1' 

# Sample model values 
Y = yFunc(X)
varY = var(Y)

# Variance gradient 
Gx = vgSA(yFunc,X,Y)
print(Gx)

# Change ux^2 by p*100%
p = 0.1
ux1 = ux*(1+p)^0.5
X1 = matrix(runif(nMC,min=-ux1*3^0.5,max=ux1*3^0.5),ncol=1) 
colnames(X1) = 'x1' 
varY1 = var(yFunc(X1))

cat('Variance change\n')
cat('MC: ',(varY1-varY)/varY,'\n')
cat('VG: ', p*Gx[1])