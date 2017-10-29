# First_Step
# First_Step
##### Exonic Content of transcripts
# exon length/transcript length
getwd()
setwd("/home/boris/Documents/First_Step/Concatenation/First_BedtoolTrial")

Transcript_Exons_Oncogene <- read.table("/home/boris/Documents/First_Step/Concatenation/First_BedtoolTrial/Transcripts_Interest_Oncogene20.fasta", sep="\t")
Transcript_Exons_Oncogene[,3] <- nchar(as.character(Transcript_Exons_Oncogene[,2]))
colnames(Transcript_Exons_Oncogene) <- c("Transcript_ID", "Sequence", "Exon_Length")


Transcript_Exons_TumorSuppressor <- read.table("/home/boris/Documents/First_Step/Concatenation/First_BedtoolTrial/Transcripts_Interest_TumorSuppressor20.fasta", sep="\t")
Transcript_Exons_TumorSuppressor[,3] <- nchar(as.character(Transcript_Exons_TumorSuppressor[,2]))
colnames(Transcript_Exons_TumorSuppressor) <- c("Transcript_ID", "Sequence", "Exon_Length")


#From the Transcripts.version ID we redownloaded the transcripts start and end from biomart
Transcript_Coord_Oncogene <- read.table("/home/boris/Documents/First_Step/Concatenation/First_BedtoolTrial/List_Transcript_Version_ToGetTranscriptLength/Transcripts_Version_Length_Oncogene.csv", sep="\t", header=T)
####The order in Transcripts_Version_Length_Oncogene.csv compared to the fasta is not the same so we will have to reorder
temp1 <- data.frame(Transcript_Coord_Oncogene[,1],Transcript_Coord_Oncogene[,3] - Transcript_Coord_Oncogene[,2])
colnames(temp1) <- c("Transcript_ID", "Transcript_Length")
#vector with the non version genes in the right order
ordered_og_List <- sapply(strsplit(as.character(Transcript_Exons_Oncogene$Transcript_ID),"\\."), `[`, 1)
temp1 <- temp1[match(ordered_og_List, temp1$Transcript_ID),]

Transcript_Exons_Oncogene[,4] <- temp1[,2]
colnames(Transcript_Exons_Oncogene) <- c("Transcript_ID", "Sequence", "Exon_Length", "Transcript_Length")


####The order in Transcripts_Version_Length_TumorSuppressor.csv compared to the fasta is not the same so we will have to reorder
Transcript_Coord_TumorSuppressor <- read.table("/home/boris/Documents/First_Step/Concatenation/First_BedtoolTrial/List_Transcript_Version_ToGetTranscriptLength/Transcripts_Version_Length_TumorSuppressor.csv", sep="\t", header=T)
temp2 <- data.frame(Transcript_Coord_TumorSuppressor[,1],Transcript_Coord_TumorSuppressor[,3] - Transcript_Coord_TumorSuppressor[,2])
colnames(temp2) <- c("Transcript_ID", "Transcript_Length")
#vector with the non version genes in the right order
ordered_TS_List <- sapply(strsplit(as.character(Transcript_Exons_TumorSuppressor$Transcript_ID),"\\."), `[`, 1)
temp2 <- temp2[match(ordered_TS_List, temp2$Transcript_ID),]

Transcript_Exons_TumorSuppressor[,4] <- temp2[,2]
colnames(Transcript_Exons_TumorSuppressor) <- c("Transcript_ID", "Sequence", "Exon_Length", "Transcript_Length")

#adding exonic content
Transcript_Exons_Oncogene[,5] <- Transcript_Exons_Oncogene[,3] / Transcript_Exons_Oncogene[,4]
colnames(Transcript_Exons_Oncogene) <- c("Transcript_ID", "Sequence", "Exon_Length", "Transcript_Length", "Exonic_Content")

Transcript_Exons_TumorSuppressor[,5] <- Transcript_Exons_TumorSuppressor[,3] / Transcript_Exons_TumorSuppressor[,4]
colnames(Transcript_Exons_TumorSuppressor) <- c("Transcript_ID", "Sequence", "Exon_Length", "Transcript_Length", "Exonic_Content")

#So since the transcript length is higher for oncogenes but that the number of exons per
#transcript and the total exon length is smaller in oncogenes, we expect a lower exonic
#content of the oncogene transcripts
#H0: The exonic content of the transcripts is the same between oncogene and tumor suppressor
#lncRNAs
#H1: The exonic content of the transcripts is different between oncogene and tumor suppressor
#lncRNAs

#calculation of mean, median, var, sd
mean(Transcript_Exons_Oncogene[,5])#0.1854489
median(Transcript_Exons_Oncogene[,5])#0.05758593
var(Transcript_Exons_Oncogene[,5])#0.0729505
sd(Transcript_Exons_Oncogene[,5])#0.2700935

mean(Transcript_Exons_TumorSuppressor[,5])#0.2039044
median(Transcript_Exons_TumorSuppressor[,5])#0.1178667
var(Transcript_Exons_TumorSuppressor[,5])#0.06495905
sd(Transcript_Exons_TumorSuppressor[,5])#0.2548707

#Checking the distribution of the data
par(mfrow=c(1,3))
qqnorm(log(Transcript_Exons_Oncogene[,5]))
qqline(log(Transcript_Exons_Oncogene[,5]))
boxplot(log(Transcript_Exons_Oncogene[,5]))
hist(log(Transcript_Exons_Oncogene[,5]))
dev.off()
#Normal when log10 transformation of the data
#The distribution of the data of the oncogene group is not normal

par(mfrow=c(1,3))
qqnorm(log(Transcript_Exons_TumorSuppressor[,5]))
qqline(log(Transcript_Exons_TumorSuppressor[,5]))
boxplot(log(Transcript_Exons_TumorSuppressor[,5]))
hist(log(Transcript_Exons_TumorSuppressor[,5]))
dev.off()
#The distribution of the data of the tumor suppressor group is normal after log10 transformation

#since the log transformation of the data is normal, will apply a t test
t.test(log(Transcript_Exons_Oncogene[,5]), log(Transcript_Exons_TumorSuppressor[,5]))
#wilcox.test(Transcript_Exons_Oncogene[,5], Transcript_Exons_TumorSuppressor[,5]) #gives the same
#significant difference, as expected the oncogene exonic content is lower since by
#definition the exonic content is: total exonic length/transcript length
#and so since the exonic length was small for oncogene and their transcript length bigger



setwd("/home/boris/Documents/First_Step/R")

codes3 <- c(as.vector(Transcript_Exons_Oncogene$Transcript_ID), as.vector(Transcript_Exons_TumorSuppressor$Transcript_ID))
df44 <- data.frame(rep(c("Oncogenes","Tumor_Suppressors"), c(454,223)), codes3, c(Transcript_Exons_Oncogene$Exonic_Content, Transcript_Exons_TumorSuppressor$Exonic_Content))
colnames(df44) <- c("File_Number", "Transcript_ID", "Exonic_Content")
#boxplot
library(ggplot2)
options(scipen=999)
pdf("Exonic_Content_of_Transcript.pdf")
ggplot(df44, aes(x=File_Number, y=Exonic_Content, fill=File_Number)) +
  labs(x="LncRNA group", y="Exonic Content") +
  ggtitle("Comparison of Oncogene and Tumor \nsuppressor LncRNAs Exonic Content of Transcripts") +
  scale_fill_discrete(name="LncRNA Type") +
  theme(legend.title = element_text(colour="black", size=15)) +#, face="bold"
  theme(legend.text = element_text(colour="black", size=13 )) +#face="bold" to get the names in bold
  theme(axis.text = element_text(size = 13))+#axis textt size
  theme(axis.title = element_text(size = 16))+#axis title size
  theme(plot.title = element_text(size = 18))+
  geom_boxplot() +
  ylim(0,1)
#geom_text(data = aggregate(Gene_Length ~ File_Number, df42, mean), aes(label = Gene_Length, y = Gene_Length + 0.08)) 
dev.off()

#same but removing some outliers

pdf("Exonic_Content_of_Transcript_Zoomed.pdf")
ggplot(df44, aes(x=File_Number, y=Exonic_Content, fill=File_Number)) +
  labs(x="LncRNA group", y="Exonic Content") +
  ggtitle("Comparison of Oncogene and Tumor \nsuppressor LncRNAs Exonic Content of Transcripts") +
  scale_fill_discrete(name="LncRNA Type") +
  theme(legend.title = element_text(colour="black", size=15)) +#, face="bold"
  theme(legend.text = element_text(colour="black", size=13 )) +#face="bold" to get the names in bold
  theme(axis.text = element_text(size = 13))+#axis textt size
  theme(axis.title = element_text(size = 16))+#axis title size
  theme(plot.title = element_text(size = 18))+
  geom_boxplot() +
  ylim(0,0.5)
#geom_text(data = aggregate(Gene_Length ~ File_Number, df42, mean), aes(label = Gene_Length, y = Gene_Length + 0.08)) 
dev.off()


#Violin plot with zoom
pdf("Exonic_Content_of_Transcript_Violin_Zoomed.pdf")
ggplot(df44, aes(x=File_Number, y=Exonic_Content, fill=File_Number)) +
  labs(x="LncRNA group", y="Exonic Content") +
  ggtitle("Comparison of Oncogene and Tumor \nsuppressor LncRNAs Exonic Content of Transcripts") +
  scale_fill_discrete(name="LncRNA Type") +
  theme(legend.title = element_text(colour="black", size=15)) +#, face="bold"
  theme(legend.text = element_text(colour="black", size=13 )) +#face="bold" to get the names in bold
  theme(axis.text = element_text(size = 13))+#axis textt size
  theme(axis.title = element_text(size = 16))+#axis title size
  theme(plot.title = element_text(size = 18))+
  ylim(0,1)+
  geom_violin() +
  geom_boxplot(width=0.1)
dev.off()

#to create a color gradient from green to blue and separate in two plots
pdf("Exonic_Content_of_Transcript_Histogram.pdf")
ggplot(df44, aes(as.vector(Exonic_Content), fill = ..x..) ) +#change x by y if the y is defined to make a vertical gradient
  geom_histogram( binwidth = 0.01) +
  scale_fill_gradient("Exonic Content",low = "blue", high = "red") +
  labs(x="Exonic Content", y="Occurence") +
  ggtitle("Comparison of Oncogene and Tumor \nsuppressor LncRNAs Exonic Content of Transcripts") +
  theme(legend.title = element_text(colour="black", size=15)) +#, face="bold"
  theme(legend.text = element_text(colour="black", size=13 )) +
  facet_grid(File_Number ~ .)
dev.off()

#let's order the values of the barplot => create a new data frame with ordered values
df44_ordered <- df44[order(df44[,3], df44[,1], decreasing=TRUE),]
df44_ordered <- data.frame(df44_ordered$File_Number, df44_ordered$Transcript_ID, df44_ordered$Exonic_Content, stringsAsFactors = FALSE)
#Here we had to set stringAsFactor as FALSE to be able to keep the order
#make the column with the values an ordered factor to keep the order in the plot
colnames(df44_ordered) <- c("File_Number", "Transcript_ID", "Exonic_Content")
df44_ordered$Transcript_ID <- factor(df44_ordered$Transcript_ID, levels = df44_ordered$Transcript_ID)


pdf("Exonic_Content_of_Transcript_Barplot.pdf")
ggplot(data=df44_ordered, aes(x=Transcript_ID, y=Exonic_Content, fill=File_Number)) +
  geom_bar(stat="identity") +
  labs(x="LncRNA group", y="Exonic Content") +
  ggtitle("Comparison of Oncogene and Tumor \nsuppressor LncRNAs Exonic Content of Transcripts") +
  theme(axis.text.x=element_blank()) +#this is to remove the names
  theme(legend.title = element_text(color="black", size=15)) +#, face="bold"
  theme(legend.text = element_text(color="black", size=13 )) +
  scale_fill_manual("LncRNAs Transcripts",values = c("firebrick3", "turquoise3")) +
  geom_hline(aes(yintercept=mean(Transcript_Exons_Oncogene[,5])), color="firebrick3", lty=1) +
  geom_hline(aes(yintercept=mean(Transcript_Exons_TumorSuppressor[,5])), color="turquoise3", lty=1)+
  theme(axis.text = element_text(size = 13))+#axis textt size
  theme(axis.title = element_text(size = 15))+#axis title size
  theme(plot.title = element_text(size = 18))
dev.off()

