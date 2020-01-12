# Create file containing all possible triples (i.e. parent 1 --> hybrid <-- parent 2) for HyDe

# Use:
# Simply need a vector of populations 
# populations <- c("pop1","pop2","pop3","pop4","pop5")
# output <- create_hyde_triples(populations)
# write.table(output, "./triples.txt", sep="\t", row.names = F, col.names = F, quote=F)


create_hyde_triples<-function(populations) {
  triples<-expand.grid(populations,populations,populations)
  dups<-vector()
  for(i in 1:nrow(triples)){
    if(any(duplicated(as.character(triples[i,]))==TRUE)){
      dups<-c(dups,i)
      }
    else{}
  }
  triples<-triples[-dups,]
  unique_triples<-matrix(ncol=3)
  for(i in populations){
    tmp<-subset(triples,triples$Var2==i)
    tmp2<-as.matrix(tmp[!duplicated(t(apply(tmp,1,sort))),])
    unique_triples<-rbind(unique_triples,tmp2)
  }
  unique_triples<-unique_triples[-1,]
  return(unique_triples)
}
