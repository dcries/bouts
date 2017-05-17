library(foreign)

temp <- read.xport("\\\\my.files.iastate.edu\\Users\\dcries\\Downloads\\DEMO_C.xpt")
temp2 <- read.xport("\\\\my.files.iastate.edu\\Users\\dcries\\Downloads\\DEMO_D.xpt")

smk <- read.xport("\\\\my.files.iastate.edu\\Users\\dcries\\Downloads\\SMQ_C.xpt")
smk2 <- read.xport("\\\\my.files.iastate.edu\\Users\\dcries\\Downloads\\SMQ_D.xpt")

temp <- temp[,c("SEQN","DMDEDUC2")]
temp2 <- temp2[,c("SEQN","DMDEDUC2")]
smk <- smk[,c("SEQN","SMQ040")]
smk2 <- smk2[,c("SEQN","SMQ040")]


temp3 <- rbind(temp,temp2)
smk3 <- rbind(smk,smk2)
names(temp3)[1] <- "ID"
names(smk3)[1] <- "ID"

nhanes <- read.csv("\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\data\\NHANES_accel.csv")

nhanes2 <- inner_join(temp3,nhanes)
nhanes3 <- inner_join(smk3,nhanes2)

names(nhanes3) <- tolower(names(nhanes3))
names(nhanes3)[2:3] <- c("smoke","education")
names(nhanes3)[7] <- c("gender")

nhanes3$smoke[nhanes3$smoke==3] <- 0
nhanes3$smoke[nhanes3$smoke==2 | nhanes3$smoke==1] <- 1
nhanes3$smoke <- nhanes3$smoke + 1 # to make same format as bouts

nhanes3$gender[nhanes3$gender==2] <- 0
nhanes3$gender[nhanes3$gender==1] <- 1
nhanes3$gender <- nhanes3$gender + 1 # to make same format as bouts

nhanes3$black <- 0
nhanes3$black[nhanes3$race==4] <- 1
nhanes3$hispanic <- 0
nhanes3$hispanic[nhanes3$race==1|nhanes3$race==2] <- 1

#remove race == 5
nhanes3 <- nhanes3[nhanes3$race != 5,]
#remove education > 5
nhanes3 <- nhanes3[nhanes3$education <= 5,]
nhanes3$education[nhanes3$education <= 4] <- 0
nhanes3$education[nhanes3$education ==5] <- 1

#ensure all covariates measured
comp <- complete.cases(nhanes3[,c("age","gender","bmi","education","black","hispanic","smoke","modvigbouts","modvigmetsb","stratum","psu","smplwt")])

nhanes4 <- nhanes3[comp,]

#ensure replicates
a=nhanes4 %>% group_by(id) %>% summarise(length(age))
nhanes5 <- nhanes4[nhanes4$id %in% a$id[which((a$`length(age)`)>1)],]

#nhanes demographics data only for matching
match <- nhanes5[which(!duplicated(nhanes5$id)),c("age","gender","bmi","education","black","hispanic","smoke","id")]

write.csv(nhanes5,file="\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\data\\NHANES_final.csv",row.names=FALSE)
write.csv(match,file="\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\data\\NHANES_match.csv",row.names=FALSE)

