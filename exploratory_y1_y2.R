library(ggplot2)
library(dplyr)

setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("bouts.csv")

y1 = bouts %>% group_by(id,rep) %>% summarise(nbouts = max(bout) ,
                                              age=age[1],
                                              gender=gender[1],
                                              bmi=bmi[1],
                                              education=education[1],
                                              employed=employed[1],
                                              income=income[1],
                                              black=black[1],
                                              hispanic=hispanic[1],
                                              smoke=smoke[1],
                                              occupation=occupation[1],
                                              married=married[1])

y2p = bouts %>% group_by(id,rep,bout) %>% summarise(totalp = sum(mets)-30,
                                                    age=age[1],
                                                    gender=gender[1],
                                                    bmi=bmi[1],
                                                    education=education[1],
                                                    employed=employed[1],
                                                    income=income[1],
                                                    black=black[1],
                                                    hispanic=hispanic[1],
                                                    smoke=smoke[1],
                                                    occupation=occupation[1],
                                                    married=married[1])
y2p$totalp[y2p$totalp==-30] <- 0
y2 <- y2p %>% group_by(id,rep) %>% summarise(total = sum(totalp),
                                             age=age[1],
                                             gender=gender[1],
                                             bmi=bmi[1],
                                             education=education[1],
                                             employed=employed[1],
                                             income=income[1],
                                             black=black[1],
                                             hispanic=hispanic[1],
                                             smoke=smoke[1],
                                             occupation=occupation[1],
                                             married=married[1])
#y2 <- y2[y2$total!=0,]



ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + theme_bw()
ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + theme_bw()
#by gender
ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~gender) + theme_bw()
ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~gender) + theme_bw()
#by education
ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~education)  + theme_bw()
ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~education) + theme_bw()
#smoke
ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~smoke) + theme_bw()
ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~smoke) + theme_bw()
#employed
ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~employed) + theme_bw()
ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~employed) + theme_bw()
#married
ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~married) + theme_bw()
ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~married) + theme_bw()
#black
ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~black) + theme_bw()
ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~black) + theme_bw()
#hispanic
ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~hispanic) + theme_bw()
ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~hispanic) + theme_bw()


#all
plot(ecdf(y2p$totalp))
#rep 1
plot(ecdf(y2p$totalp[y2p$rep==1]))
plot(ecdf(y2p$totalp[y2p$rep==2]))
# diff
aa <- y2p %>% group_by(id,rep) %>% summarise(s=sum(totalp))
aaa <- aa %>% group_by(id) %>% filter(length((id))==2) %>% summarise(diff=s[2]-s[1])
plot(ecdf(aaa$diff))


#y1 vs y2
plot(y1$nbouts,y2$total)
plot(y1$nbouts,y1$age)
plot(y1$nbouts,y1$bmi)
plot(y2$total,y2$bmi)
plot(y2$total,y2$age)
