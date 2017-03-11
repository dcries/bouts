#plots for dissertation
library(ggplot2)
library(dplyr)
library(reshape)
library(rjags)
library(scales)
library(rcompanion)
library(reshape)

setwd("C:\\Users\\dcries\\github\\bouts\\data")
boutsraw <- read.csv("bouts.csv")
bouts <- read.csv("finalbouts.csv")

#plot 1
#showing the data, y1=nbouts and ~y2 total excess minutes
#size 4x6
q1 <- qplot(data=bouts,x=nbouts,geom="bar")  + theme_bw() + xlab("Number of Bouts in a 24 Hour Period") + theme(axis.title=element_text(size=9))
q2 <- qplot(data=bouts,x=total) + theme_bw() + xlab("Total MET-mins in MVPA in a 24 Hour Period")+ theme(axis.title=element_text(size=9))
grid.arrange(q1,q2,nrow=1)

#plot 1
#size 5 x 6
#by gender
bouts$gender2 <- factor(bouts$gender,labels=c("Female","Male"),levels=c(1,2))
q1f <- ggplot(data=subset(bouts,!is.na(gender2)),aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~gender2) + theme_bw() + xlab("Number of Bouts") + ylab("Proportion")
q1m <- ggplot(data=subset(bouts,!is.na(gender2)),aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~gender2) + theme_bw() + xlab("Total MET-mins in MVPA in a 24 Hour Period") + ylab("Density")
grid.arrange(q1f,q1m,nrow=2)

#plot 2
# plot y1 vs y2 boxplots
#size 5x6
bouts$avgtotalexcess <- 0 
bouts$avgtotalexcess[bouts$nbouts != 0] <- bouts$totalexcess[bouts$nbouts != 0]/bouts$nbouts[bouts$nbouts != 0]

q3 <- qplot(data=bouts,x=(nbouts),y=totalexcess,geom="boxplot",group=as.factor(nbouts)) +
  xlim(c(0.5,13.5)) + theme_bw() + xlab("Number of Bouts in 24 Hours") + 
  ylab("Total Excess MET-mins in 24 Hours")+ theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))

q4 <- qplot(data=bouts,x=nbouts,y=avgtotalexcess,geom="boxplot",group=nbouts)+
  xlim(c(0.5,13.5)) + theme_bw() + xlab("Number of Bouts in 24 Hours") + 
  ylab("Average Excess MET-mins by Number of Bouts in 24 Hours")+ theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))

grid.arrange(q3,q4,nrow=1)

#plot 3
# 5x6
# y1 and y2 faceted by weekend
bouts$Weekend2 <- factor(bouts$Weekend,labels=c("Weekday","Weekend"),levels=c(0,1))
q5 <- ggplot(data=subset(bouts,!is.na(Weekend2)),aes(x=nbouts,y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
  geom_bar() + facet_wrap(~Weekend2) + theme_bw() + ylab("proportion") + xlab("Number of Bouts") #+ ggtitle("faceted by weekend(1) or weekday(0)")

q6 <- ggplot(data=subset(bouts,!is.na(Weekend2)),aes(x=totalexcess,y=..density..)) + 
  geom_histogram() + facet_wrap(~Weekend2) + theme_bw() + xlab("Total Excess Minutes")

grid.arrange(q5,q6,nrow=2)



#compare std dev of obs grouped by gender, BMI, age, etc
#size 3.5 x 6.5
bouts2rep <- bouts %>% group_by(id) %>% filter(length(id)==2)
trues <- bouts2rep %>% group_by(id) %>% filter(length(id)==2) %>% summarise(s=sd(nbouts),
                                                                             s2=sd(totalexcess))
full <- cbind(trues,bouts2rep[bouts2rep$rep==1,c("bmi","age","gender")])
agegroup <- quantile(full$age,probs=c(0.25,0.5,0.75),na.rm=TRUE)
full$a <- 0
full$a[full$age > agegroup[1]] <- 1
full$a[full$age > agegroup[2]] <- 2
full$a[full$age > agegroup[3]] <- 3
bmigroup <- quantile(full$bmi,probs=c(0.25,0.5,0.75),na.rm=TRUE)
full$b <- 0
full$b[full$bmi > bmigroup[1]] <- 1
full$b[full$bmi > bmigroup[2]] <- 2
full$b[full$bmi > bmigroup[3]] <- 3

full$a2 <- factor(full$a,labels=c(paste0("<",agegroup[1]),paste0(agegroup[1],"-",agegroup[2]),
                                  paste0(agegroup[2],"-",agegroup[3]),paste0(">",agegroup[3])),
                  levels=0:3)
full$b2 <- factor(full$a,labels=c(paste0("<",bmigroup[1]),paste0(bmigroup[1],"-",bmigroup[2]),
                                  paste0(bmigroup[2],"-",bmigroup[3]),paste0(">",bmigroup[3])),
                  levels=0:3)

full$gender2 <- factor(full$gender,labels=c("Male","Female"),levels=c(2,1))

q7a <- qplot(data=full,y=(s+.0),x=(a2),group=(a2),geom="boxplot") + theme_bw() + ylab("Std Dev of Nbouts within Person") + xlab("Age Groups") + theme(axis.text.x=element_text(size=8,angle=90))
q7b <- qplot(data=full,y=(s+.0),x=(b2),group=(b2),geom="boxplot") + theme_bw()  + ylab("Std Dev of Nbouts within Person") + xlab("BMI Groups")+ theme(axis.text.x=element_text(size=8,angle=90))
q7c <- qplot(data=subset(full,!is.na(gender2)),y=(s+.0),x=gender2,group=(gender2),geom="boxplot") + theme_bw()  + ylab("Std Dev of Nbouts within Person") + xlab("Gender")+ theme(axis.text.x=element_text(size=8,angle=90))

grid.arrange(q7a,q7b,q7c,nrow=1)

#same as previous plot except for totalexcess
#size 5x8
q8a <- qplot(data=full,y=(s2+.0),x=(a2),group=(a2),geom="boxplot") + theme_bw() + ylab("Std Dev of Total Excess MET-mins within Person") + xlab("Age Groups")+ theme(axis.text.x=element_text(size=8,angle=90),axis.title.y=element_text(size=8))
q8b <- qplot(data=full,y=(s2+.0),x=(b2),group=(b2),geom="boxplot") + theme_bw()  + ylab("Std Dev of Total Excess MET-mins within Person") + xlab("BMI Groups")+ theme(axis.text.x=element_text(size=8,angle=90),axis.title.y=element_text(size=8))
q8c <- qplot(data=subset(full,!is.na(gender2)),y=(s2+.0),x=gender2,group=(gender2),geom="boxplot") + theme_bw()  + ylab("Std Dev of Total Excess MET-mins within Person") + xlab("Gender")+ theme(axis.text.x=element_text(size=8,angle=90),axis.title.y=element_text(size=8))
grid.arrange(q8a,q8b,q8c,nrow=1)



#2 way table plot
meas2 <- unique(bouts$id)[which(table(bouts$id)==2)]
bouts2rep <- bouts[bouts$id %in% meas2,]
bouts2rep$nbouts[bouts2rep$nbouts >= 11] <- 11 #cap nbouts at 20

heat <- bouts2rep %>% group_by(id) %>% summarise(nbouts1=nbouts[1],
                                                 nbouts2=nbouts[2],
                                                 out=paste0(nbouts[1]," ",nbouts[2]),
                                                 diff=nbouts[1]-nbouts[2])

count <- rep(0,121)
for(i in 0:10){
  for(j in 0:10){
    count[i*11+j+1] <- sum(heat$nbouts1==i & heat$nbouts2==j)
  }
}


heatdf <- data.frame(nbouts1=rep(0:10,each=11),nbouts2=rep(0:10,11),count=count)
qn <- quantile(heatdf$count,probs=c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))
qn01 <- rescale(c(qn, range(heatdf$count))) 

ggplot(data=heatdf,aes(x=as.factor(nbouts1),y=as.factor(nbouts2))) + geom_tile(aes(fill=count)) + 
  geom_text(aes(label=count)) + #+ scale_fill_gradient(low = "white", high = "red") 
  scale_fill_gradientn (
    colours = colorRampPalette(c( "white", "red"))(120),
    values = c(0,qn01[2:7],1)) + xlab("Number of Bouts Trial 1") + ylab("Number of Bouts Trial 2")

