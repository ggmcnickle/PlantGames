#Gordon G. McNickle et al 2014
#University of Illinois at CHicago, Bio Sci
#Wilfrid Laurier University, Biology
#gmcnickle@wlu.ca

rm(list=ls(all=TRUE))
library(scatterplot3d)
library(rootSolve)
graphics.off()

###############################################
#Parameters
###############################################
## Rn    = Root resource, N concentration uM
## Rc    = Shoot resource, CO2 concentration ppm
## u[i]  = Leaf, stem or root production strategy
## c1rn  = root cost in N units, player 1 
## c1rc	 = root cost in C units, player 1
## c1wn	 = wood cost in N units, player 1
## c1wc  = wood cost in C units, player 1
## c1ln	 = leaf cost in N units, player 1
## c1lc  = leaf cost in C units, player 1
## a1l	 = Encounter rate between leaf and C, player 1
## a1r	 = Encounter rate between root and N, player 1
## zw	 = size asymmetry constant for height competition
##		Larger numbers of z make the size-asymmetry more severe
## zl	 = size asymmetry constant for leaf competition
## alph1 = represents C:N ratio. Or plant need for C relative to N, player 1 
## beta1 = 1-alph1, player 1	
################################################

plant.game<-function(u) {
		with (as.list(params),	{
	#NB: u[1] is root, u[2] is wood, u[3] is leaf for plant 1
	#and u[4] is root, u[5] is wood, u[6] is leaf for plant 2

		#BELOWGROUND GAME		
	    	r   <- a1r*u[1]+a2r*u[4]	#Total roots
		Hn  <- Rn*(1-exp(-r))	#harvest of N by both plants
		f1  <- a1r*u[1]/(r) 	#P1 share of harvest
		f2  <- a2r*u[4]/(r)		#P2 share of harvest
		
		#DERIVATIVES	
		dHn1 <- a1r*Rn*exp(-r) 	#Derivative of Hn with respect to u (either)
		dHn2 <- a2r*Rn*exp(-r) 	#Derivative of Hn with respect to u (either)
		df1 <- a1r/r-(a1r^2)*u[1]/(r^2)	#Derivative of f1 with respect to u1
		df2 <- a2r/r-(a2r^2)*u[4]/(r^2)	#Derivative of f2 with respect to u2

		#ABOVEGROUND GAME
	    	l   <- a1l*u[3]+a2l*u[6]		#Total leaves
		lz  <- ((a1l*u[3])^zl)*(u[2]^zw) + ((a2l*u[6])^zl)*(u[5]^zw)
		Hc  <- Rc*(1-exp(-l))		#harvest of L by both plants
		z1  <- (((a1l*u[3])^zl)*(u[2]^zw))/lz		#P1 share of harvest
		z2  <- (((a2l*u[6])^zl)*(u[5]^zw))/lz		#P2 share of harvest
		
		#DERIVATIVES	
		dHc1   <- a1l*Rc*exp(-l) 	#Derivative of Hc with respect to u (either)
		dHc2  <- a2l*Rc*exp(-l) 	#Derivative of Hc with respect to u (either)

		dz1.l <- (zl*a1l*(u[2]^zw)*(u[3]^(zl-1))*lz - a1l*(u[2]^zw)*(u[3]^zl)*zl*a1l*(u[2]^zw)*(u[3]^(zl-1)))/(lz^2)	#Derivative of z1 with respect to u3
		dz2.l <- (zl*a1l*(u[5]^zw)*(u[6]^(zl-1))*lz - a1l*(u[5]^zw)*(u[6]^zl)*zl*a1l*(u[5]^zw)*(u[6]^(zl-1)))/(lz^2)	#Derivative of z2 with respect to u6
		dz1.w <- (zw*a1l*(u[3]^zl)*(u[2]^(zw-1))*lz - a1l*(u[3]^zl)*(u[2]^zw)*zw*a1l*(u[3]^zl)*(u[2]^(zw-1)))/(lz^2)	#Derivative of z1 with respect to u2
		dz2.w <- (zw*a1l*(u[6]^zl)*(u[5]^(zw-1))*lz - a1l*(u[6]^zl)*(u[5]^zw)*zw*a1l*(u[6]^zl)*(u[5]^(zw-1)))/(lz^2)	#Derivative of z2 with respect to u5

		#PROFIT FUNCTIONS
		P1n <- f1*Hn - c1rn*u[1] - c1wn*u[2] - c1ln*u[3]  #Net profit N P1
		P2n <- f2*Hn - c2rn*u[4] - c2wn*u[5] - c2ln*u[6]  #Net profit N P2
		P1c <- z1*Hc - c1rc*u[1] - c1wc*u[2] - c1lc*u[3]  #Net profit c P1
		P2c <- z2*Hc - c2rc*u[4] - c2wc*u[5] - c2lc*u[6]  #Net profit c P2

		dP1n.dur <- df1*Hn + f1*dHn1 - c1rn
		dP1n.duw <- -c1wn
		dP1n.dul <- -c1ln

		dP2n.dur <- df2*Hn + f2*dHn2 - c2rn
		dP2n.duw <- -c2wn
		dP2n.dul <- -c2ln

		dP1c.dur <- -c1rc
		dP1c.duw <- dz1.w*Hc -c1wc
		dP1c.dul <- dz1.l*Hc + z1*dHc1 -c1lc

		dP2c.dur <- -c2rc
		dP2c.duw <- dz2.w*Hc -c2wc
		dP2c.dul <- dz2.l*Hc + z2*dHc2 -c2lc
				
	#derivatives of whole plant G-function where, G=(Ps^alph)*(Pr^beta)
		dG1.dur <- alph1*(P1c^(alph1-1)) * (P1n^beta1)*dP1c.dur +
			  beta1*(P1c^alph1) * (P1n^(beta1-1))*dP1n.dur
		dG1.duw <- alph1*(P1c^(alph1-1)) * (P1n^beta1)*dP1c.duw +
			  beta1*(P1c^alph1) * (P1n^(beta1-1))*dP1n.duw
		dG1.dul <- alph1*(P1c^(alph1-1)) * (P1n^beta1)*dP1c.dul +
			  beta1*(P1c^alph1) * (P1n^(beta1-1))*dP1n.dul

		dG2.dur <- alph2*(P2c^(alph2-1)) * (P2n^beta2)*dP2c.dur +
			  beta2*(P2c^alph2) * (P2n^(beta2-1))*dP2n.dur
		dG2.duw <- alph2*(P2c^(alph2-1)) * (P2n^beta2)*dP2c.duw +
			  beta2*(P2c^alph2) * (P2n^(beta2-1))*dP2n.duw
		dG2.dul <- alph2*(P2c^(alph2-1)) * (P2n^beta2)*dP2c.dul +
			  beta2*(P2c^alph2) * (P2n^(beta2-1))*dP2n.dul


		return(c(dG1.dur = dG1.dur, dG1.duw = dG1.duw, dG1.dul = dG1.dul, 
			dG2.dur = dG2.dur, dG2.duw = dG2.duw, dG2.dul = dG2.dul))
		})} 

#Loop to randomize traits
cc.dist1<-c(60:120)/100 	#uniform distribuion of C costs from .6 to 3
cc.dist<-c(90:250)/100 	#uniform distribuion of C costs from .6 to 3
a.dist<-c(400:1000)/1000
#THIS C RANGE MAKES RESPIRATION 50% OF GPP
cn.dist<-c(10:140)/1000 		#uniform distribuion of N costs from 0.1 to 0.3
alpha.dist<-c(8000:9900)/10000	#uniform distribuion of alphas from 0.95 to 0.999	
carb<-c(50:1500)		#Carbon availability from 1000 - 5000
nitr<-c(25:80)		#Nitrogen availability from 500 -2000
Zl.dist<-(100:150)/100
Zw.dist<-(110:500)/100

#Create empty vectors to store output
Root1<-as.numeric()	#empty root vector
Wood1<-as.numeric()	#empty wood vector
Leaf1<-as.numeric()	#empty leaf vector
Root2<-as.numeric()	#empty root vector
Wood2<-as.numeric()	#empty wood vector
Leaf2<-as.numeric()	#empty leaf vector
prec<-as.numeric()	#empty precision vector
cost<-as.numeric()
n.cost<-as.numeric()
alpha.vec1<-as.numeric()
alpha.vec2<-as.numeric()
N.vec<-as.numeric()


imax<-35000				#Max number of randomizations
i<-0

#LOOP
while (i<imax) {
		a<-sample(a.dist, 4, replace=TRUE)
		cr<-sample(cc.dist1, 6, replace=TRUE)
		c<-sample(cc.dist, 6, replace=TRUE)
		n<-sample(cn.dist, 6, replace=TRUE)
		alp<-sample(alpha.dist, 2, replace=TRUE)
		R.c<-sample(carb, 1, replace=TRUE)
		R.n<-sample(nitr, 1, replace=TRUE)
		zl<-sample(Zl.dist, 1, replace=TRUE)
		zw<-sample(Zw.dist, 1, replace=TRUE)

		params<- c(Rc=R.c, Rn=R.n, zl=zl, zw=zw, 
				a1r=a[1], a1l=a[2], c1rc=cr[1], c1wc=cr[2], c1lc=c[3], c1rn=n[1], c1wn=n[1], c1ln=n[1], alph1=alp[1], beta1=(1-alp[1]), 
				a2r=a[3], a2l=a[4], c2rc=cr[4], c2wc=cr[5], c2lc=c[6], c2rn=n[4], c2wn=n[4], c2ln=n[4], alph2=alp[2], beta2=(1-alp[2]))
		solution<-multiroot(f = plant.game, start = c(.1,.1,.1,.1,.1,.1), maxiter=1000, positive = TRUE)
		ESS<-solution$root
		
		#if (sum(ESS)=0, 
			
		prec <- solution$estim.precis
		
		Root1<-c(Root1,if (!is.nan(prec)) {if (prec < 1e-6) {ESS[1]} else {0}} else {0})
		Wood1<-c(Wood1,if (!is.nan(prec)) {if (prec < 1e-6) {ESS[2]} else {0}} else {0})
		Leaf1<-c(Leaf1,if (!is.nan(prec)) {if (prec < 1e-6) {ESS[3]} else {0}} else {0})
		Root2<-c(Root2,if (!is.nan(prec)) {if (prec < 1e-6) {ESS[4]} else {0}} else {0})
		Wood2<-c(Wood2,if (!is.nan(prec)) {if (prec < 1e-6) {ESS[5]} else {0}} else {0})
		Leaf2<-c(Leaf2,if (!is.nan(prec)) {if (prec < 1e-6) {ESS[6]} else {0}} else {0})
		
		y<-cr[1]*ESS[1]+cr[2]*ESS[2]+c[3]*ESS[3]+cr[4]*ESS[4]+cr[5]*ESS[5]+c[6]*ESS[6]
		y.n<-n[1]*ESS[1]+n[2]*ESS[2]+n[3]*ESS[3]+n[4]*ESS[4]+n[5]*ESS[5]+n[6]*ESS[6]
		cost<-c(cost, if (min(ESS)>0) {y} else {0})
		n.cost<-c(n.cost, if (min(ESS)>0) {y.n} else {0})
		alpha.vec1<-c(alpha.vec1, alp[1])
		alpha.vec2<-c(alpha.vec2, alp[2])
		N.vec<-c(N.vec, R.n)
		i<-i+1
		print(i)
		flush.console()  # force the output 
		}

Root<-alpha.vec1*Root1+alpha.vec2*Root2
Leaf<-alpha.vec1*Leaf1+alpha.vec2*Leaf2
Wood<-alpha.vec1*Wood1+alpha.vec2*Wood2

Rootx<-Root1+Root2
Leafx<-Leaf1+Leaf2
Woodx<-Wood1+Wood2
Totalx<-Rootx+Leafx+Woodx

Total<-Root+Leaf+Wood
fRoot<-(Root)/Total
fWood<-(Wood)/Total
fLeaf<-(Leaf)/Total

GPP<-Total+cost

out<-data.frame(Root, Wood, Leaf, Total, fRoot, fWood, fLeaf, cost, alpha.vec1, alpha.vec2, N.vec, GPP)
write.csv(out, "C:/Users/gmcnickle/Desktop/out.3T_TEMP.csv") 

dev.new()
par(mfrow=c(3,2))
plot((Leaf)~Total, xlim=c(0,1400), ylim=c(0,700))
plot(fWood~fRoot, xlim=c(0,1), ylim=c(0,1))
plot((Wood)~Total, xlim=c(0,1400), ylim=c(0,700))
plot(fLeaf~fRoot, xlim=c(0,1), ylim=c(0,1))
plot((Root)~Total, xlim=c(0,1400), ylim=c(0,700))
plot(fWood~fLeaf, xlim=c(0,1), ylim=c(0,1))

dev.new()
Total[Total<1]<-NA
Totalx[Totalx<1]<-NA
hist(Total)

CUE.reg<-lm(Total~GPP)
b<-round(as.numeric(CUE.reg$coefficients[1]),2)
b.t<-if(b<0) {b} else {paste("+",b)}
m<-round(as.numeric(CUE.reg$coefficients[2]),2)
text<-paste("NPP=",m,"*GPP",b.t)
 
dev.new()
par(mfrow=c(1,2), mar=c(5,5,2,2), cex.axis=1.4, cex.lab=1.9)
plot(GPP, Total, xlab="GPP", ylab="NPP", xlim=c(0,4000), ylim=c(0,2000))
mtext("(a)", side=3, adj=0.05, line=-1.3, cex=1.5)
mtext(text, side=3, adj=0.25, line=-4, cex=1.25)
abline(a = b, b = m, lty = 3, lwd=2)
#boxplot(na.omit(Total/(cost+Total)), ylab="Carbon use efficiency", ylim=c(0,.85))
#mtext("(b)", side=3, adj=0.05, line=-1.3)
boxplot(na.omit((Totalx)/N.vec), ylab="NPP per N available", ylim=c(0,60))
mtext("(b)", side=3, adj=0.05, line=-1.3, cex=1.5)




