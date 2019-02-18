# belowground turnover rate (y^-1)
BGTR<-1
# rootshoot ratio
RS<-2
# organic self-packing density
k1<-.085
# inorganic self-packing density
k2<-1.99
# stable fraction of organic matter
kr<-.1
# number of tides in a year
f<-704
# capture rate of suspended sediment
q<-2.8
# maximum biomass
ymax<-2000
# depth of minimun biomass (rel to MHW)
minD<--30
# mean high water, make permutations
MHW<-seq(5, 120, by=5)
# mean low water
MLW<--MHW
# k is a counter
k<-0
# permutations of depth
D<-seq(-30,130,by= 2)
# define the length of the output arrays
dsdt<-numeric(length(MHW)*length(D))
dodt<-numeric(length(MHW)*length(D))
dzdt<-numeric(length(MHW)*length(D))
NEC<-numeric(length(MHW)*length(D))
Dseq<-numeric(length(MHW)*length(D))
MHWseq<-numeric(length(MHW)*length(D))
# iterate through the permutations of MHW
for (i in 1:length(MHW)){
  maxD<-MHW[i]+10
maxE<-MHW[i]-minD
Dopt<-(maxD+minD)/2
# compute the biomass coefficients
a <- -((-minD * ymax - maxD * ymax) / ((minD - Dopt) * (-maxD + Dopt)))
b <- -(ymax / ((minD - Dopt) * (-maxD + Dopt)))
c <- (minD * maxD * ymax) / ((minD - Dopt) * (maxD - Dopt))
# iterate throught the permutions of depth
for (j in 1:length(D)){
# compute the biomass
 Bs<-( a * D[j] + b * (D[j]) ^ 2 + c) 
 if (Bs<0) Bs<-0 
# compute the marsh elevation for each permutation of MHW and D
Z<- MHW[i]  -D[j]
# FIT is the fractional inuntation time
FIT <-(MHW[i] - Z) / (MHW[i] - MLW[i])
if (FIT>1) FIT<-1
if (FIT<0) FIT<-0
# the product of FIT and q should not be greater than 1
# i.e. capture of sediment in any single tide must be <=1
qstar<- q * FIT
if (q > 1) qstar <- 1
if (q < 0) qstar<- 0
k<-k+1
Ds<-D[j]
if(D[j]<0 ) Ds<-0
Dseq[k]=D[j]
MHWseq[k]=MHW[i]
# compute the accretion rates from mineral and organic inputs
dsdt[k] <- (0.5 * Ds * 0.000021 * f * qstar) / k2
dodt[k] <- kr * RS * BGTR * 0.0001 * Bs / k1
# compute the total vertical accretion rate
dzdt[k] <- dsdt[k] + dodt[k]
# NEC is Normalized Elevation Capital
#'NEC =  (Z-(MSL-10cm))/((MHW+30cm)-(MSL-10cm))
NEC[k] = (Z + 10) / (MHW[i] + 30 + 10)
# the output is limited to 0<NEC<1
if (NEC[k]<0) NEC[k]<-NA
 }

}
# split the vectors into pieces of the length of D for plotting
NECsplit<-split(NEC, ceiling(seq_along(NEC)/length(D)))
Dsplit<-split(Dseq, ceiling(seq_along(Dseq)/length(D)))                
dzsplit<-split(dzdt, ceiling(seq_along(dzdt)/length(D)))
dosplit<-split(dodt,ceiling(seq_along(dodt)/length(D)))
dssplit<-split(dsdt,ceiling(seq_along(dsdt)/length(D)))
par(new=F)
# plot total mineral accretion for permulations of D against NEC
plot(NECsplit[[1]],dssplit[[1]],"l",xlim=c(0,1),ylim=c(0,.8),xlab='NEC',ylab=' Accretion (cm/y)', col=2)
  par(new=T)
for (i in 2:length(MHW)){
  plot(NECsplit[[i]],dssplit[[i]],"l",xlim=c(0,1),ylim=c(0,.8),xlab='',ylab='',axes=F, col=2)
  par(new=T)}
# plot total accretion for permulations of D against NEC
for (i in 1:length(MHW)){
plot(NECsplit[[i]],dzsplit[[i]],"l",xlim=c(0,1),ylim=c(0,.8),xlab='',ylab='',axes=F, col=1)
par(new=T)
}
# plot total organic accretion for permulations of D against NEC
for (i in 1:length(MHW)){
  plot(NECsplit[[i]],dosplit[[i]],"l",xlim=c(0,1),ylim=c(0,.8),xlab='',ylab='',axes=F, col=3)
  par(new=T)
}
# stip the NA cases 
dataout<-data.frame(MHWseq,NEC,Dseq,dsdt,dodt,dzdt) 
finalcases<-dataout[complete.cases(dataout),]
# write the good solutions to a CSV file
write.csv(finalcases,file="MEMout.csv",col.names=TRUE,sep=" ")

