setwd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/StorageHyL")
files<-dir(pattern="outerrHybUlt")

out<-c()
for(i in 1:length(files)){
Fi<-files[i]
logi1<-grep(x=readLines(Fi)[3],pattern="Dot figure generated in file")
if(length(logi1)==0){a<-c()}else {
a<-substring(gsub(pattern=".txt", replacement="",Fi), 13)}
out<-c(out,a)
}

setwd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne")
destinoHy<-"/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/UltraMetric/Ult_JuliaForHybLam/"
destinoJu<-"/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/UltraMetric/Ult_RJuliaNetworks/"
destinoR<-"/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/UltraMetric/Ult_RNetworks/"

Hy1<-"/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/JuliaForHybLam/"
Ju1<-"/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/RJuliaNetworks/"
R1<-"/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/RNetworks/"

hyFi<-paste("JulHybrLamb",out,sep="")
rFi<-paste("RNetwork",out,sep="")
juFi<-paste("RJuliaNet",out,".txt",sep="")



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