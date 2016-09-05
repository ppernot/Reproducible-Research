# Modèle
fExpr = function(x1) 1 + exp(-x1^2)*cos(x1*10)

x1 = 0.5; names(x1)='x1'

# Pdf
x1.pdf = 'norm'

# Range of variation of uX
uMin=0.01; uMax=4; npt = 100

# MC runs
nMC = 10000 

