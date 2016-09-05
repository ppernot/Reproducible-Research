dir.create(path='figs', 
           showWarnings = FALSE)

#### 
model= 'ExpUnif'

load(file=paste0('results/',model,'.Rda'))
source(paste0('models/',model,'.R'))

png(file=paste0('figs/fig',model,'0.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
xx=seq(-10,10,by=0.001)
yy=exp(-abs(xx))
zz=xx*0; zz[abs(xx)<=2]=0.25
plot(xx,yy,xlab='X',ylab='Y, PDF(X)',
     type='l',lty=1,lwd=4,col=4,cex=1,
     xlim=c(-5,5),ylim=c(0,1.1),yaxs='i')
lines(xx,zz,lty=1,col=2,lwd=4)
legend(1.5,1,c('PDF(X)',expression(Y == e^{-abs(X)})),
       col=c(2,4),lwd=4,bty='n',cex=2.5)
abline(v=0,lty=2)
box(lwd=2)
dev.off()

png(file=paste0('figs/fig',model,'1.png'),width=1200,height=900)
par(mar=c(5,5,3,1),cex=1.5,cex.lab=2,cex.axis=2,cex.main=2)
sel = seq(1,nrow(res),by=2)
xx=res[sel,1]; yy=res[sel,2]; zz=res[sel,4]
plot(xx,yy,xlab=expression(u[X]),ylab='Arbitrary units',
     main=expression(Y == e^{-abs(X)}),
     type='p',pch=15,col=4,cex=2,
     ylim=c(-0.5,1),yaxs='i',
     xlim=c(0,6),xaxs='i')
points(xx,zz,pch=17,col=2,cex=2)
points(xx,res[sel,5],pch=19,col=1,cex=2)
legend(2.5,0.8,c(expression(paste(u[Y],' (MC)')),
                 expression(paste(G[X],' (MC)')),
                 expression(paste(IQ[95],' (MC)')),
                 expression(paste(u[Y],' (Ana)')),
                 expression(paste(G[X],' (Ana)')) ),
       col=c(4,2,1,4,2),lty=c(0,0,0,2,2),lwd=4,
       pch=c(15,17,19,-1,-1),pt.cex=2,bty='n',cex=2,
       ncol=2, title='X ~ Unif', title.adj=0)
abline(v=1.91,lwd=2,lty=2) # Position du maximum
abline(h=0,lty=2,lwd=2)
# Analytical for Unif distrib
a=xx*3^0.5
mY = (1-exp(-a))/a
varY=(1-exp(-2*a))/(2*a)-mY^2
lines(xx,varY^0.5,col=4,lwd=4,lty=2)
gX=(0.5*(2*a+1)*exp(-2*a) -1/2 - 2*mY *((a+1)*exp(-a)-1)) / (2*a*varY)
lines(xx,gX,col=2,lwd=4,lty=2)
box(lwd=2)
dev.off()

png(file=paste0('figs/fig',model,'2.png'),width=1200,height=1800)
par(mfrow=c(3,1),mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2,cex.main=2)
for ( a in c(1,3.3,8) ) {
  x1.u=a/3^0.5
  S=rgumlib::gumS1(nrun=50000,adapt=FALSE,silent=TRUE,
                   fExpr=fExpr,x.mu=x1,x.u=x1.u,x.pdf=x1.pdf)
  vg = vgSA(fExpr=fExpr,x.mu=x1,x.u=x1.u,X=S$X, Y=S$Y)$vg
  hist(S$Y,col='gold',xlim=c(0,1),xlab='Y', freq=FALSE, 
       main=paste0('uX = ', signif(x1.u,2),'; ',
                   'uY = ', signif(sd(S$Y),2),'; ',
                   'Gx = ', signif(vg,2))
  )
  
}
dev.off()

nMC=50000 
npt=1000 
uMin=0.01; uMax=4.01
dx=(uMax-uMin)/(npt-1)
breaks=seq(0,1,by=0.01)
res=matrix(0,ncol=length(breaks)-1,nrow=npt)
for ( i in 1:npt ) {
  x1.u = uMin +(i-1)*dx
  S=rgumlib::gumS1(nrun=nMC,adapt=FALSE,silent=TRUE,
                   fExpr=fExpr,x.mu=x1,x.u=x1.u,x.pdf=x1.pdf)
  res[i,] = hist(S$Y, breaks=breaks,plot=FALSE)$counts 
}


png(file=paste0('figs/fig',model,'_test.png'),width=2400,height=1200)
par(mfrow=c(1,2),mar=c(6,6,1.5,0),cex=1.5,cex.lab=2,cex.axis=2,cex.main=2)
cols=fields::two.colors(512,start="white",middle="yellow",end="red")
x=seq(uMin,uMax,length.out=nrow(res))
y=0.5*(breaks[-length(breaks)]+breaks[-1])
show=c(0.1,1.0,1.9,3.0)
image(y,x,t(res),col=cols,ylab=expression(u[X]),xlab='Y',
      yaxs='i',ylim=c(0,4),xaxt='n')
axis(1,lwd=3,at=seq(0.1,0.9,by=0.1))
contour(y,x,t(res),col='black',add=TRUE,nlevels=25,drawlabels=FALSE)
abline(h=show,lty=2,lwd=3,col='black')
# abline(h=1.9 ,lty=2,lwd=3,col='green')
box(lwd=3)

par(mar=c(6,1,1.5,6))
plot(y,y,type='n',,xaxs='i',yaxt='n',xlab='Y',yaxs='i',ylim=c(0,4),bty='n')
axis(1,lwd=3,at=seq(0.1,0.9,by=0.1))
grid(ny=0,lty=2,lwd=3,col='brown')
abline(v=c(0.005,0.995),lty=2,lwd=3,col='brown')
abline(h=show,lty=2,lwd=3,col='black')
# abline(h=1.9 ,lty=2,lwd=3,col='green')

rmax=max(res)
for(u in show){
  i = which(x>u)[1]
  segments(y,u+0.01,y,u+0.01+res[i,]/rmax*uMax,lwd=12,
           col='orchid',lend=2)  
}

a=show*3^0.5
mY = (1-exp(-a))/a
varY=(1-exp(-2*a))/(2*a)-mY^2
uY=varY^0.5
gX=(0.5*(2*a+1)*exp(-2*a) -1/2 - 2*mY *((a+1)*exp(-a)-1)) / (2*a*varY)
lab = paste0('uY = ',signif(uY,2))
text(0.7,y=show+0.4,labels=lab,cex=2)
lab = paste0('Gx = ',signif(gX,1))
text(0.7,y=show+0.6,labels=lab,cex=2)

# box(lwd=3)
# mtext(show,side=4,at=show,adj=1,las=1)
dev.off()

###########################################################################
model= 'ExpNorm'

load(file=paste0('results/',model,'.Rda'))
source(paste0('models/',model,'.R'))

png(file=paste0('figs/fig',model,'0.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
xx=seq(-10,10,by=0.001)
yy=exp(-abs(xx))
zz=dnorm(xx,mean=0,sd=2/3)
plot(xx,yy,xlab='X',ylab='Y, g(X)',
     type='l',lty=1,lwd=4,col=4,cex=1,
     xlim=c(-5,5),ylim=c(0,1.1),yaxs='i')
lines(xx,zz,lty=1,col=2,lwd=4)
legend(1.5,1,c('PDF(X)',expression(Y == e^{-abs(X)})),
       col=c(2,4),lwd=4,bty='n',cex=2.5)
abline(v=0,lty=2)
box(lwd=2)
dev.off()


png(file=paste0('figs/fig',model,'1.png'),width=1200,height=900)
par(mar=c(5,5,3,1),cex=1.5,cex.lab=2,cex.axis=2,cex.main=2)
sel = seq(1,nrow(res),by=2)
xx=res[sel,1]; yy=res[sel,2]; zz=res[sel,4]
plot(xx,yy,xlab=expression(u[X]),
     ylab='Arbitrary units',
     main = expression(Y == e^{-abs(X)}),
     type='p',pch=15,col=4,cex=2,
     ylim=c(-0.5,1),yaxs='i',
     xlim=c(0,6),xaxs='i')
points(xx,zz,pch=17,col=2,cex=2)
points(xx,res[sel,5],pch=19,col=1,cex=2)
legend(2.5,0.8,c(expression(paste(u[Y],' (MC)')),
                 expression(paste(G[X],' (MC)')),
                 expression(paste(IQ[95],' (MC)')) ),
       col=c(4,2,1),lty=c(0,0,0),lwd=4,
       pch=c(15,17,19),pt.cex=2,bty='n',cex=2,
       ncol=1, title='X ~ Norm', title.adj=0)
abline(v=1.9,lwd=2,lty=2) # Position du maximum
abline(h=0,lty=2,lwd=2)

# # Analytical for Unif distrib
# a=xx*3^0.5
# mY = (1-exp(-a))/a
# varY=(1-exp(-2*a))/(2*a)-mY^2
# lines(xx,varY^0.5,col=4,lwd=4)
# gX=(0.5*(2*a+1)*exp(-2*a) -1/2 - 2*mY *((a+1)*exp(-a)-1)) / (2*a*varY)
# lines(xx,gX,col=2,lwd=4)
box(lwd=2)
dev.off()

###########################################################################
model= 'SinUnif'

load(file=paste0('results/',model,'.Rda'))
source(paste0('model/',model,'.R'))

png(file=paste0('figs/fig',model,'0.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
xx=seq(-2*pi,2*pi,by=0.001)
yy=fExpr(xx)
zz=xx*0; zz[abs(xx)<=2]=0.25
plot(xx,yy,xlab='X',ylab='Y, PDF(X)',
     type='l',lty=1,lwd=4,col=4,cex=1,
     xlim=c(-2*pi,2*pi),xaxs='i',
     ylim=c(-1.1,1.1),yaxs='i')
lines(xx,zz,lty=1,col=2,lwd=4)
legend(-4,1,c('PDF(X)',expression(Y == sin(X))),
       col=c(2,4),lwd=4,bty='n',cex=2.5)
abline(v=0,lty=2,lwd=2); abline(h=0,lty=2,lwd=2)
box(lwd=2)
dev.off()

png(file=paste0('figs/fig',model,'1.png'),width=1200,height=900)
par(mar=c(6,6,3,1.5),cex=1.5,cex.lab=2,cex.axis=2)
sel = seq(1,nrow(res),by=2)
xx=res[sel,1]; yy=res[sel,2]; zz=res[sel,4]
plot(xx,yy,xlab=expression(u[X]),ylab='Arbitrary units',
     type='p',pch=15,col=4,cex=2,
     ylim=c(-0.7,2.2),yaxs='i',
     xlim=c(0,uMax),xaxs='i')
points(xx,zz,pch=17,col=2,cex=2)
points(xx,res[sel,5],pch=19,col=1,cex=2)
abline(v=pi/3^0.5*(c(0.5,0.75,1,1.25,1.5,1.75)),
       lwd=2,lty=2,col='brown') # Position du maximum
abline(h=0,lty=2,lwd=2)
abline(a=0,b=1,lty=1,lwd=4,col=4)
legend(1.8,1.8,c(expression(paste(u[Y],' (MC)')),
               expression(paste(u[Y],' (GUM)')),
               expression(G[X]),
               expression(IQ[95]),
               expression(paste(u[Y],' (Ana)')),
               expression(paste(G[X],' (Ana)')) ),
       col=c(4,4,2,1,4,2),lty=c(0,1,0,0,2,2),lwd=4,
       pch=c(15,-1,17,19,-1,-1),pt.cex=2,bty='o',cex=2,
       ncol=2, title='X ~ Unif', title.adj=0, bg='white',box.col='white')

# Analytical
x=c(); y=c(); z=c(); i=0
for (ux in seq(0,5,by=0.01)) {
  a=ux*3^0.5
  varYT=0.5-0.25*sin(2*a)/a
  GxT = 0.25*(sin(2*a)/(2*a) -cos(2*a))/varYT
  i=i+1; x[i]=ux; y[i]=varYT^0.5; z[i]=GxT
}
lines(x,y,col='blue',lwd=4,lty=2)
lines(x,z,col='red',lwd=4,lty=2)

box(lwd=2)
mtext(c('a=',expression(pi/2),expression(3*pi/4),
        expression(pi),expression(5*pi/4),
        expression(3*pi/2),expression(7*pi/4)),
      side=3,padj=0.5,line=1.5,
      at=c(0.7,pi/3^0.5*c(0.5,0.75,1,1.25,1.5,1.75)),cex=3)
dev.off()


png(file=paste0('figs/fig',model,'2.png'),width=1200,height=1800)
par(mfrow=c(3,1),mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2,cex.main=2)
for ( a in c(0.5,2.25,3^0.5*uMax) ) {
  x1.u= a/3^0.5
  S=rgumlib::gumS1(nrun=50000,adapt=FALSE,silent=TRUE,
                   fExpr=fExpr,x.mu=x1,x.u=x1.u,x.pdf=x1.pdf)
  vg = vgSA(fExpr=fExpr,x.mu=x1,x.u=x1.u,X=S$X, Y=S$Y)$vg
  hist(S$Y,col='gold',xlim=c(-1,1),xlab='Y', freq=FALSE, 
       main=paste0('uX = ', signif(x1.u,2),'; ',
                   'uY = ', signif(sd(S$Y),2),'; ',
                   'Gx = ', signif(vg,2))
  )
}
dev.off()

# png(file=paste0('figs/fig',model,'3.png'),width=1200,height=1800)
nMC=50000# 50000
npt=200
uMin=0.01; uMax=4.01
dx=(uMax-uMin)/(npt-1)
breaks=seq(-1,1,by=0.02)
res1=matrix(0,ncol=length(breaks)-1,nrow=npt)
for ( i in 1:npt ) {
  x1.u = uMin +(i-1)*dx
  S=rgumlib::gumS1(nrun=nMC,adapt=FALSE,silent=TRUE,
                   fExpr=fExpr,x.mu=x1,x.u=x1.u,x.pdf=x1.pdf)
  res1[i,] = hist(S$Y, breaks=breaks,plot=FALSE)$counts 
}
# png(file=paste0('figs/fig',model,'3.png'),width=1200,height=1200)
# par(mfrow=c(1,1),mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2,cex.main=2)
# cols=fields::two.colors(512,start="white",middle="yellow",end="red")
# x=seq(uMin,uMax,length.out=nrow(res))
# y=0.5*(breaks[-length(breaks)]+breaks[-1])
# image(x,y,res,col=cols,xlab=expression(u[X]),ylab='Y')
# contour(x,y,res,col='black',add=TRUE,nlevels=21,drawlabels=FALSE)
# abline(v=1.3,lwd=2,col='green')
# box(lwd=2)
# dev.off()

png(file=paste0('figs/fig',model,'_test.png'),width=2400,height=1200)
par(mfrow=c(1,2),mar=c(6,6,1.5,0),cex=1.5,cex.lab=2,cex.axis=2,cex.main=2)
cols=fields::two.colors(512,start="white",middle="yellow",end="red")
x=seq(uMin,uMax,length.out=nrow(res1))
y=0.5*(breaks[-length(breaks)]+breaks[-1])
show=c(0.2,1.4,2.25,3.1)
image(y,x,t(res1),col=cols,ylab=expression(u[X]),xlab='Y',
      yaxs='i',ylim=c(0,4),xaxt='n')
axis(1,lwd=3,at=seq(-0.9,0.9,by=0.2))
contour(y,x,t(res1),col='black',add=TRUE,nlevels=55,drawlabels=FALSE)
abline(h=show,lty=2,lwd=3,col='black')
# abline(h=1.9 ,lty=2,lwd=3,col='green')
box(lwd=3)

par(mar=c(6,1,1.5,6))
plot(y,y,type='n',,xaxs='i',yaxt='n',xlab='Y',yaxs='i',ylim=c(0,4),bty='n')
axis(1,lwd=3,at=seq(-0.9,0.9,by=0.2))
grid(ny=0,lty=2,lwd=3,col='brown')
abline(v=c(-0.992,0.992),lty=2,lwd=3,col='brown')
abline(h=show,lty=2,lwd=3,col='black')
# abline(h=1.9 ,lty=2,lwd=3,col='green')

rmax=max(res1)
uY = gX = c()
j=0
for(u in show){
  j=j+1
  i = which(x>u)[1]
  segments(y,u+0.01,y,u+0.01+res1[i,]/rmax*uMax,lwd=12,
           col='orchid',lend=2) 
  S=rgumlib::gumS1(nrun=nMC,adapt=FALSE,silent=TRUE,
                   fExpr=fExpr,x.mu=x1,x.u=u,x.pdf=x1.pdf)
  uY[j] = sd(S$Y)
  gX[j]=rgumlib::vgSA(fExpr=fExpr,x.mu=x1,x.u=u,S$X,S$Y,silent=TRUE,budgetTable=FALSE)$vg
}

lab = paste0('uY = ',signif(uY,2))
text(0.25,y=show+0.4,labels=lab,cex=2)
lab = paste0('Gx = ',signif(gX,1))
text(0.25,y=show+0.6,labels=lab,cex=2)

# box(lwd=3)
# mtext(show,side=4,at=show,adj=1,las=1)
dev.off()

###########################################################################
model= 'SinNorm'

load(file=paste0('results/',model,'.Rda'))
source(paste0('models/',model,'.R'))

png(file=paste0('figs/fig',model,'1.png'),width=1200,height=900)
par(mar=c(6,6,3,1.5),cex=1.5,cex.lab=2,cex.axis=2)
sel = seq(1,nrow(res),by=2)
xx=res[sel,1]; yy=res[sel,2]; zz=res[sel,4]
plot(xx,yy,xlab=expression(u[X]),ylab='Arbitrary units',
     type='p',pch=15,col=4,cex=2,
     ylim=c(-0.7,2.2),yaxs='i',
     xlim=c(0,uMax),xaxs='i')
points(xx,zz,pch=17,col=2,cex=2)
points(xx,res[sel,5],pch=19,col=1,cex=2)
abline(v=pi/3^0.5*(c(0.5,0.75,1,1.25,1.5,1.75)),lwd=2,lty=2, col='brown') # Position du maximum
abline(h=0,lty=2,lwd=2)
abline(a=0,b=1,lty=1,lwd=4,col=4)
legend(2,1.8,c(expression(paste(u[Y],' (MC)')),
               expression(paste(u[Y],' (GUM)')),
               expression(G[X]),
               expression(IQ[95])),
       col=c(4,4,2,1),lty=c(0,1,0,0),lwd=4,
       pch=c(15,-1,17,19),pt.cex=2,bty='o',cex=2,
       ncol=2, title='X ~ Norm', title.adj=0, bg='white',box.col='white')
box(lwd=2)
mtext(c('a=',expression(pi/2),expression(3*pi/4),
        expression(pi),expression(5*pi/4),
        expression(3*pi/2),expression(7*pi/4)),
      side=3,padj=0.5,line=1.5,
      at=c(0.7,pi/3^0.5*c(0.5,0.75,1,1.25,1.5,1.75)),cex=3)
dev.off()

###########################################################################
model= 'TanhUnif'

load(file=paste0('results/',model,'.Rda'))
source(paste0('model/',model,'.R'))

png(file=paste0('figs/fig',model,'0.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
xx=seq(-5,5,by=0.001)
yy=fExpr(xx)
zz=dunif(xx,-1,1)
plot(xx,yy,xlab='X',ylab='Y, PDF(X)',
     type='l',lty=1,lwd=4,col=4,cex=1,
     ylim=c(-1.1,1.1),yaxs='i')
lines(xx,zz,lty=1,col=2,lwd=4)
legend('topleft',c('PDF(X)',expression(Y == tanh(X))),
       col=c(2,4),lwd=4,bty='n',cex=2.5)
abline(v=0,lty=2,lwd=2); abline(h=0,lty=2,lwd=2)
box(lwd=2)
dev.off()

png(file=paste0('figs/fig',model,'1.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
sel = seq(1,nrow(res),by=2)
xx=res[sel,1]; yy=res[sel,2]; zz=res[sel,4]
plot(xx,yy,xlab=expression(u[X]),ylab='Arbitrary units',
     type='p',pch=15,col=4,cex=2,
     ylim=c(0,2.1),yaxs='i',
     xlim=c(0,uMax),xaxs='i')
points(xx,zz,pch=17,col=2,cex=2)
points(xx,res[sel,5],pch=19,col=1,cex=2)
legend(4,1.8,c(expression(paste(u[Y],' (MC)')),
               expression(paste(u[Y],' (GUM)')),
               expression(G[X]),
               expression(IQ[95])),
       col=c(4,4,2,1),lty=c(0,1,0,0),lwd=4,
       pch=c(15,-1,17,19),pt.cex=2,bty='n',cex=2,
       ncol=1)
j=attr(eval(deriv(expression(tanh(x)),'x'),list(x=x1)),'gradient')
abline(a=0,b=j,lty=1,lwd=4,col=4)
# abline(v=1.88,lwd=2,lty=2) # Position du maximum
abline(h=c(0,1,2),lty=2,lwd=2)
box(lwd=2)
dev.off()

###########################################################################
model= 'TanhUnif2'

load(file=paste0('results/',model,'.Rda'))
source(paste0('model/',model,'.R'))

png(file=paste0('figs/fig',model,'0.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
xx=seq(-6,6,by=0.001)
yy=fExpr(xx)
zz=dunif(xx,0.5,2.5)
plot(xx,yy,xlab='X',ylab='Y, PDF(X)',
     type='l',lty=1,lwd=4,col=4,cex=1,
     ylim=c(-1.1,1.1),yaxs='i')
lines(xx,zz,lty=1,col=2,lwd=4)
legend('topleft',c('PDF(X)',expression(Y == tanh(X))),
       col=c(2,4),lwd=4,bty='n',cex=2.5)
abline(v=0,lty=2,lwd=2); abline(h=0,lty=2,lwd=2)
box(lwd=2)
dev.off()

png(file=paste0('figs/fig',model,'1.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
xx=res[,1]; yy=res[,2]; zz=res[,4]
plot(xx,yy,xlab=expression(u[X]),ylab='Arbitrary units',
     type='p',pch=19,col=4,cex=1,
     ylim=c(0,2.1),yaxs='i',
     xlim=c(0,uMax),xaxs='i')
points(xx,zz,pch=19,col=2,cex=1)
points(xx,res[,5],pch=19,col=1,cex=1.5)
legend(4,1.8,c(expression(paste(u[Y],' (MC)')),
               expression(paste(u[Y],' (GUM)')),
               expression(G[X]),
               expression(IQ[95])),
       col=c(4,4,2,1),lty=c(0,1,0,0),lwd=4,
       pch=c(19,-1,19,19),pt.cex=1.5,bty='n',cex=2,
       ncol=1)
j=attr(eval(deriv(expression(tanh(x)),'x'),list(x=x1)),'gradient')
abline(a=0,b=j,lty=1,lwd=4,col=4)
# abline(v=1.88,lwd=2,lty=2) # Position du maximum
abline(h=c(0,1,2),lty=2,lwd=2)
box(lwd=2)
dev.off()
###########################################################################
model= 'DampCosNorm'

load(file=paste0('results/',model,'.Rda'))
source(paste0('model/',model,'.R'))

png(file=paste0('figs/fig',model,'0.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
xx=seq(-4,6,by=0.001)
yy=fExpr(xx)
zz=dnorm(xx,mean=x1,sd=2/3)
plot(xx,yy,xlab='X',ylab='Y, g(X)',
     type='l',lty=1,lwd=4,col=4,cex=1,
     ylim=c(0,2.1),yaxs='i')
lines(xx,zz,lty=1,col=2,lwd=4)
legend(1,2,c('g(X)','Y'),
       col=c(2,4),lwd=4,bty='n',cex=2.5,
       title=expression(Y == 1 + exp(-x^2)*cos(10*x)))
abline(v=0,lty=2)
box(lwd=2)
dev.off()

png(file=paste0('figs/fig',model,'1.png'),width=1200,height=900)
par(mar=c(6,6,1.5,1.5),cex=1.5,cex.lab=2,cex.axis=2)
xx=res[,1]; yy=res[,2]; zz=res[,4]
plot(xx,yy,xlab=expression(u[X]),ylab='Arbitrary units',
     type='p',pch=19,col=4,cex=1,
     ylim=c(-0.7,1),yaxs='i',
     xlim=c(0,2),xaxs='i')
points(xx,zz,pch=19,col=2,cex=1)
legend(6,0.9,c(expression(u[Y]),expression(G[X])),
       col=c(4,2),lwd=4,pch=19,pt.cex=1,bty='n',cex=2.5,
       title=expression(Y == e^{-abs(X)}))
# abline(v=1.88,lwd=2,lty=2) # Position du maximum
abline(h=0,lty=2,lwd=2)
box(lwd=2)
dev.off()

