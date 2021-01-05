# Create file containing all possible pairwise comparisons
# Use:
# Simply need a vector of populations 
# populations <- c("pop1","pop2","pop3","pop4","pop5")
# output <- create_pairwise_comparisons(populations)

create_pairwise_comparisons<-function(populations) {
  doubles<-expand.grid(populations,populations)
  dups<-vector()
  for(i in 1:nrow(doubles)){
    if(any(duplicated(as.character(doubles[i,]))==TRUE)){
      dups<-c(dups,i)
    }
    else{}
  }
  doubles<-doubles[-dups,]
  unique_doubles<-doubles[!duplicated(t(apply(doubles,1,sort))),]
  return(unique_doubles)
}
