

install.packages("igraph")
library("igraph")


kmerize <- function(dna_sequences, k) {
  kmers <- c()
  for (seq in dna_sequences) {
    for (i in 1:(nchar(seq) - k + 1)) {
      kmers <- c(kmers, substr(seq, i, i + k - 1))
    }
  }
  kmers
}

kmerGraph <- function(kmers) {
  edges <- c()
  for (i in 1:length(kmers)) {
    for (j in 1:length(kmers)) {
      if (substr(kmers[i], 2, nchar(kmers[i])) == substr(kmers[j], 1, nchar(kmers[j]) - 1)) {
        edges <- c(edges, c(kmers[i], kmers[j]))
      }
    }
  }
  graph_from_edgelist(matrix(edges, ncol = 2, byrow = TRUE))
}

# Function to check if a graph is Eulerian
isEulerian <- function(graph) {
  degree_seq <- degree(graph)
  all(degree_seq %% 2 == 0)
}

getAssembly <- function(graph) {
  if (!isEulerian(graph)) return(NULL)
  eulerian_path <- eulerian(graph)
  sequence <- ""
  for (i in 1:length(eulerian_path$vertex)) {
    sequence <- paste0(sequence, substr(eulerian_path$vertex[i], 1, 1))
  }
  sequence <- paste0(sequence, substr(eulerian_path$vertex[length(eulerian_path$vertex)], 2, nchar(eulerian_path$vertex[length(eulerian_path$vertex)])))
  sequence
}

readWriter <- function(dna_sequence, readSize, numReads) {
  reads <- c()
  for (i in 1:numReads) {
    start_pos <- sample(1:(nchar(dna_sequence) - readSize + 1), 1)
    reads <- c(reads, substr(dna_sequence, start_pos, start_pos + readSize - 1))
  }
  reads
}

readWronger <- function(reads, errorRate) {
  set.seed(42)
  total_bases <- sum(nchar(reads))
  num_errors <- round(total_bases * errorRate)
  error_positions <- sample(1:total_bases, num_errors)
  for (pos in error_positions) {
    read_idx <- findInterval(pos, cumsum(nchar(reads)))
    base_pos <- pos - sum(nchar(reads)[1:(read_idx - 1)])
    original_base <- substr(reads[read_idx], base_pos, base_pos)
    possible_bases <- setdiff(c("A", "T", "C", "G"), original_base)
    new_base <- sample(possible_bases, 1)
    substr(reads[read_idx], base_pos, base_pos) <- new_base
  }
  reads
}

# Test Run
dna_sequence <- "ATTCAATTCAG"
k <- 5
kmers <- kmerize(dna_sequence, k)
graph <- kmerGraph(kmers)
assembled_sequence <- getAssembly(graph)
print(assembled_sequence)
