#20 august, gorkemer, some scripts for data wrangling


a<- unique(d$participant_ID)
a_oddDeleted <- a[-which(a==17952)]
a_oddDeleted <- a_oddDeleted[-which(a_oddDeleted==19101)]
a_oddDeleted <- a_oddDeleted[-which(a_oddDeleted==19198)]

b<- unique(data$subNum)
b_added <- append(b, c(20000))

#sort
a_oddDeleted <- sort(a_oddDeleted)
b_added <- sort(b_added)

mat <- cbind(a_oddDeleted, b_added)