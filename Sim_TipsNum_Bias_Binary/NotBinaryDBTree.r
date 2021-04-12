
#No binary tree in birth dead simulation with NetSim even usign mrca=TRUE

#In birth dead simulation you can find  two scenaries:
#tree goes extinct=0 and no extinct tips are sampled=1

#In adition to, we find that NetSim produces only not binary tree with mrca=FALSE
#mrca=FALSE produces both, binary and not binary trees
#diversitree can't stimate parametres with not binary tree so you have to remove it

#

library(NetSim) # simulate bd time-tree networks
library(TreeSim) # simulate bd time-trees
library(ape) # tools for handling phylo objects
library(diversitree) # tools for estimating lambda and mu
library("ggplot2")  # Data visualization
library("dplyr")    # Data manipulation

numbsim=500
bdNS <- NetSim::sim.bdh.age(age = 9, numbsim = numbsim, lambda = 0.9, mu = 0.5, nu = 0,hybprops = c(0.5, 0.5, 0.5), hyb.inher.fxn = NetSim::make.beta.draw(1, 1), frac = 0.9,mrca = TRUE, complete = FALSE, stochsampling = FALSE, hyb.rate.fxn = NULL, trait.model = NULL)
## remove any trees with no taxa
bdNS <- bdNS[!sapply(X = bdNS, FUN = is.null)]
#Information about tree goes extinct=0 and no extinct tips are sampled=1
bdNSNF <- bdNS[!sapply(X = bdNS, FUN = is.phylo)]#birth dead tree with NetSim with value 0 or 1
NotData<-as.data.frame(table(unlist(bdNSNF)))#Frequency of 0 or 1 tree
## remove non-phylo elements
bdNS <- bdNS[sapply(X = bdNS, FUN = is.phylo)]
##not Binary tree
RTNum<-sum(!sapply(X = bdNS, FUN = is.binary))
# remove non-binary
bdNS <- bdNS[sapply(X = bdNS, FUN = is.binary)]

#Frequency of tree with value 0,1 and binary and no binary tree.
Infor<-rbind(NotData,data.frame(Var1="Rooted",Freq=RTNum),data.frame(Var1="Complete",Freq=length(bdNS)))
Infor$prop<-Infor$Freq/numbsim



#codes to make a nice pie with ggplot2
colnames(Infor)<-c("class","n","prop")
count.data<-Infor

# Add label position
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data$Por<-count.data$prop*100
count.data

#mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

#pdf(file="rereg.pdf",width=7, height=7)
ggplot(data = count.data, 
       aes(x = 2, y = prop, fill = class))+
  geom_bar(stat = "identity")+
  coord_polar("y", start = 200) +
  geom_text(aes(y = lab.ypos, label = paste(Por,"%", sep = "")), col = "white") +
  theme_void() +
  scale_fill_brewer(palette = "Dark2")+
  xlim(.2,2.5)+ggtitle("Percentage of data to estimate parameters  \n NetSim(age=9,lambda=0.9)") 
#dev.off()