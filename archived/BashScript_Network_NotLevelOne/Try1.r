arg1<-as.numeric(commandArgs(trailingOnly=TRUE)[1])
arg2<-as.numeric(commandArgs(trailingOnly=TRUE)[2])
arg3<-as.numeric(commandArgs(trailingOnly=TRUE)[3])
arg4<-commandArgs(trailingOnly=TRUE)[4]
print(arg1)
print(arg2)
print(arg3)
res<-rnorm(arg1,arg2,arg3)
print(res)
print(boxplot(res))

eval(parse(text=as.character(arg4)))
print(numbsim1)
print(n1)
