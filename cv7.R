rm(list=ls())
setwd("V:/MPA/MPAPRG/exercise_07")

# Libraries
library(Biostrings)


## Task 1
Score <- function(s, sequences, l){
  motifs <- character(length(sequences))
  
  for (i in 1:length(sequences)) {
    motifs[i] <- substr(as.character(sequences[i]), s[i], s[i] + l - 1)
    print(motifs)
  }
  motifs <- as.matrix(motifs)
  
  score <- 0
  
  for (column in 1:l){
    occurences <- table(motifs[,column])
    score = score + as.numeric(occurences[which.max(occurences)])
    
  }
  print(score)
  return(score)
}

sequences <- DNAStringSet(c("ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"))

print(Score(c(1,1,1,1), sequences, 7))

## Task 2
NextLeaf <- function(s, t, k){
  for (i in t:1){
    if(s[i] < k){
      s[i] <- s[i] + 1
      return(s)
    }
    s[i] <- 1
    }
  return(s)
}

## Task 3
BFMotifSearch <- function(DNA, t, n, l){
  s <- rep(1, times=t)
  bestScore <- Score(s, DNA, l)
  while (TRUE){
    s <- NextLeaf(s, t, n - l + 1)
    if (Score(s, DNA, l) > bestScore){
      bestScore <- Score(s, DNA, l)
      bestMotif <- s
    }
    if (all(s == rep(1, times=t))){
      return(bestMotif)
    }
  }
}

DNA <- readDNAStringSet('seq_motif.fasta')

print(BFMotifSearch(DNA, 4, 8, 3))


## Task 4
NextVertex <- function(s, i, t, k){
  if (i < t){
    s[i + 1] <- 1
    return(s, i+1)
  }else{
    for(j in t:1){
      if(s[j] < k){
        s[j] <- s[j] + 1
        return(s, j)
      }
    }
  }
  return(s, 0)
}

## Task 5
ByPass <- function(s, i, t, k){
  for (j in i:1){
    if (s[j] < k){
      s[j] <- s[j] + 1
      return(s, j)
    }
  }
  return(s, 0)
}

## Task 6
# BBMotifSearch<- function(DNA, t, n, l){
#   s <- rep(1, t)
#   bestScore <- 0
#   i <- 1
#   while (i > 0){
#     if(i < t){
#       optimisticScore <- Score(s, i, DNA, l) + (t - i) * l
#       if (optimisticScore < bestScore){
#         (s, i) <- ByPass(s, i, t, n - l + 1)
#       }else{
#         (s, i) <- NextVertex(s, i, t, n - l + 1)
#       }
#     }else if (Score(s, t, DNA, l) > bestScore){
#       bestScore <- Score(s, t, DNA, l)
#       bestMotif <- s
#       (s, i) <- NextVertex(s, i, t, n − l + 1)
#     }else{
#       
#     }
#   }
#   return(bestMotif)
# }




















