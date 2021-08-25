arg1<-commandArgs(trailingOnly=TRUE)[1]
arg2<-commandArgs(trailingOnly=TRUE)[2]
arg3<-commandArgs(trailingOnly=TRUE)[3]
arg4<-commandArgs(trailingOnly=TRUE)[4]
arg5<-commandArgs(trailingOnly=TRUE)[5]
arg6<-commandArgs(trailingOnly=TRUE)[6]
arg7<-commandArgs(trailingOnly=TRUE)[7]
arg8<-commandArgs(trailingOnly=TRUE)[8]

arg9<-commandArgs(trailingOnly=TRUE)[9]
arg10<-commandArgs(trailingOnly=TRUE)[10]
arg11<-commandArgs(trailingOnly=TRUE)[11]
arg12<-commandArgs(trailingOnly=TRUE)[12]


#$patStoHL ->1
#$pathR ->2
#$pathRJ ->3
#$pathJHL ->4
#$PatUltR ->5
#$PatUltRJ ->6
#$PatUlJHL ->7
#$PatUltHLStor ->8

#$PatNuR 9
#$PatNuRJ 10
#$PatNuJHL 11
#$PatNuHLSt 12


setwd(arg1)
files<-dir(pattern="outerrHybUlt")


out<-c()
for(i in 1:length(files)){
Fi<-files[i]
logi1<-grep(x=readLines(Fi)[3],pattern="Dot figure generated in file")
if(length(logi1)==0){a<-c()}else {
a<-substring(gsub(pattern=".txt", replacement="",Fi), 13)}
out<-c(out,a)
}


All1<-substring(gsub(pattern=".txt", replacement="",files), 13)
NUOut<-All1[!All1%in%out]



#setwd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne")

destinoHy<-paste(arg7,"/",sep="")
destinoJu<-paste(arg6,"/",sep="")
destinoR<-paste(arg5,"/",sep="")

Hy1<-paste(arg4,"/",sep="")
Ju1<-paste(arg3,"/",sep="")
R1<-paste(arg2,"/",sep="")

hyFi<-paste("JulHybrLamb",out,sep="")
rFi<-paste("RNetwork",out,sep="")
juFi<-paste("RJuliaNet",out,".txt",sep="")



NUhyFi<-paste("JulHybrLamb",NUOut,sep="")
NUrFi<-paste("RNetwork",NUOut,sep="")
NUjuFi<-paste("RJuliaNet",NUOut,".txt",sep="")


#
HLsto<-c(paste("HybUltHybUlt",out,".pdf",sep=""),
paste("HybUltHybUlt",out,"_coal_unit",sep=""),
paste("HybUltHybUlt",out,".dot",sep=""),
paste("outerrHybUlt",out,".txt",sep=""),
paste("JulHybrLamb",out,sep=""))

targetdir <- paste(arg1,"/",HLsto,sep="")
filestocopy <- paste(arg8,"/",HLsto,sep="")

file.copy(from=targetdir, to=filestocopy,copy.mode = TRUE,recursive = FALSE)
file.remove(targetdir)

#
HLsto<-c(paste("HybUltHybUlt",NUOut,".pdf",sep=""),
paste("HybUltHybUlt",NUOut,"_coal_unit",sep=""),
paste("HybUltHybUlt",NUOut,".dot",sep=""),
paste("outerrHybUlt",NUOut,".txt",sep=""),
paste("JulHybrLamb",NUOut,sep=""))

targetdir <- paste(arg1,"/",HLsto,sep="")
filestocopy <- paste(arg12,"/",HLsto,sep="")

file.copy(from=targetdir, to=filestocopy,copy.mode = TRUE,recursive = FALSE)
file.remove(targetdir)



#


##
origindir <- Hy1
targetdir <- paste(Hy1,hyFi,sep="")
filestocopy <- paste(destinoHy,hyFi,sep="")

file.copy(from=targetdir, to=filestocopy,copy.mode = TRUE,recursive = FALSE)
file.remove(targetdir)


#####
origindir <- Ju1
targetdir <- paste(Ju1,juFi,sep="")
filestocopy <- paste(destinoJu,juFi,sep="")

file.copy(from=targetdir, to=filestocopy,copy.mode = TRUE,recursive = FALSE)
file.remove(targetdir)


#####
origindir <- R1
targetdir <- paste(R1,rFi,sep="")
filestocopy <- paste(destinoR,rFi,sep="")

file.copy(from=targetdir, to=filestocopy,copy.mode = TRUE,recursive = FALSE)
file.remove(targetdir)



#Not Ultrametric
NUhyFi<-paste("JulHybrLamb",NUOut,sep="")
NUrFi<-paste("RNetwork",NUOut,sep="")
NUjuFi<-paste("RJuliaNet",NUOut,".txt",sep="")

#$PatNuR 9
#$PatNuRJ 10
#$PatNuJHL 11

##
origindir <- Hy1
targetdir <- paste(Hy1,NUhyFi,sep="")
filestocopy <- paste(arg11,"/",NUhyFi,sep="")

file.copy(from=targetdir, to=filestocopy,copy.mode = TRUE,recursive = FALSE)
file.remove(targetdir)


#####
origindir <- Ju1
targetdir <- paste(Ju1,NUjuFi,sep="")
filestocopy <- paste(arg10,"/",NUjuFi,sep="")

file.copy(from=targetdir, to=filestocopy,copy.mode = TRUE,recursive = FALSE)
file.remove(targetdir)


#####
origindir <- R1
targetdir <- paste(R1,NUrFi,sep="")
filestocopy <- paste(arg9,"/",NUrFi,sep="")

file.copy(from=targetdir, to=filestocopy,copy.mode = TRUE,recursive = FALSE)
file.remove(targetdir)
