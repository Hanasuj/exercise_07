rm(list=ls())
setwd("E:/asdghjk/skola/VUT Brno/inzinier/4rocnik/PRG")

# Libraries
library(Biostrings)


## Task 1
Score <- function(s, DNA, l){
  alignment_matrix <- c()
  for (i in 1:length(s)){
    seq <- subseq(DNA[[i]], start = s[i], width = l)
    alignment_matrix <- c(alignment_matrix, seq)
  }
  alignment_matrix <- DNAStringSet(alignment_matrix)
  
  frequency_profile <- consensusMatrix(alignment_matrix)
  
  consensus_string <- c()
  
  score <- 0
  for (i in 1:l){
    score <- score + max(frequency_profile[,i])
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
    return(list(s, i+1))
  }else{
    for(j in t:1){
      if(s[j] < k){
        s[j] <- s[j] + 1
        return(list(s, j))
      }
    }
  }
  return(list(s, 0))
}

## Task 5
ByPass <- function(s, i, t, k){
  for (j in i:1){
    if (s[j] < k){
      s[j] <- s[j] + 1
      return(list(s, j))
    }
  }
  return(list(s, 0))
}

## Task 6
BBMotifSearch<- function(DNA, t, n, l){
  s <- rep(1, t)
  bestScore <- 0
  i <- 1
  while (i > 0){
    if(i < t){
      optimisticScore <- Score(s, DNA, l) + (t - i) * l
      if (optimisticScore < bestScore){
        bypass <- ByPass(s, i, t, n - l + 1)
        s <- bypass[[1]]
        i <- bypass[[2]]
     }else{
        next_vertex <- NextVertex(s, i, t, n - l + 1)
        s <- next_vertex[[1]]
        i <- next_vertex[[2]]
      }
    }else if (Score(s, DNA, l) > bestScore){
      bestScore <- Score(s, DNA, l)
      bestMotif <- s
      next_vertex <- NextVertex(s, i, t, n - l + 1)
      s <- next_vertex[[1]]
      i <- next_vertex[[2]]
    }else{
      
    }
  }
  return(bestMotif)
}

result <- BBMotifSearch(DNA, 4, 8, 3)

print(result)









