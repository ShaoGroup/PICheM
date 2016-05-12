# Organize the Compound-Protien Interactions (CPIs)  Data.
# Songpeng Zu /zusongpeng@gmail.com/

#-- Include packages
#library(ggplot2)
library(data.table)

#-- Information of Protein Classification from Homo sapiens
Pro2Class <- read.table("Protein2Class_Chembl.txt",sep="\t",header=FALSE,quote="",
                        col.names=c("protein","class"))

ProClassName <- read.table("ProClassName.txt",sep="\t",header=FALSE,quote="",colClasses = "character",
                           col.names=c("class","detail"))

OrganPro <- read.table("organism.txt",sep="\t",header=FALSE,quote="",colClasses=c("character","factor"),
                       col.names=c("protein","organism"))

total <- merge(Pro2Class,OrganPro,by="protein")
total <- unique(total)
# We noticed that four proteins have multiple classes from the source file. We modified them and change 
# the data of "Protein2Class_chembl.txt".
# The four proteins: CHEMBL2683 CHEMBL3521 CHEMBL3646 CHEMBL5213
# The modified results:
# CHEMBL2683 Transporter enzyme Both are right. (Do not need to change.)
# CHEMBL3521 phosphatase
# CHEMBL3646 phosphatase
# CHEMBL5213 protein kinase

pro.human <- total[which(total$organism=="Homo sapiens"),]
#pie(table(pro.human$class))
#barplot(table(pro.human$class),cex.names = 0.5,las=2)
#barplot(table(pro.human$class),cex.names = 0.5,horiz=TRUE,las=2,col="red")
# proclass.stat <- as.data.frame(table(pro.human$class))
# colnames(proclass.stat) <- c("ProteinClass","Frequency")
# bp <- ggplot(proclass.stat, aes(ProteinClass,Frequency)) + geom_bar(fill="red",stat="identity") 
# bp +theme(axis.text.y = element_text(colour = "black",size=16),
#           axis.title.y = element_blank(),axis.title.x=element_text(colour="black",
#                                                                    size=16)) +coord_flip()
#-- Information of Known Protein-Compound Interactions from Chembl
CPIs <- read.table("compounds_targets_interaction_7.txt",sep="\t",header=FALSE,quote="",
                   colClasses=c("character","character","numeric","NULL"),
                   col.names=c("compound","protein","affinity","NA"))
CPIs.class <- merge(CPIs,pro.human,by="protein")

#-- Generate the CPIsByR file.
# Only proteins from homo sapiens are included.
write.table(CPIs.class,file="CPIsByR.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

#-- Generate the CPIs data by protein family.

for(i in 1:nrow(ProClassName)){
  CPIsTemp <- CPIs.class[grep(ProClassName[i,1],CPIs.class$class),1:3]
  # Everage the binding affinities if a compound have multi records agains one protein.
  z <- as.data.table(CPIsTemp)
  keys <- colnames(z)[!grepl('affinity',colnames(z))]
  CPIsTemp.ave <- z[,list(aveaffinty=mean(affinity)),keys]
  
  filetmp <- paste("CPIs_",ProClassName[i,1],".txt")
  write.table(CPIsTemp.ave,file=filetmp,quote=FALSE,sep="\t",row.names=FALSE,col.names = TRUE)
}

