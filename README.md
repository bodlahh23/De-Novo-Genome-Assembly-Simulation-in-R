# De-Novo-Genome-Assembly-Simulation-in-R
A simplified simulation of de novo genome assembly using k-mers and Eulerian graphs in R.
This project demonstrates the core computational principles behind genome assembly using graph theory concepts similar to those used in modern assemblers.

## 📌 **Overview**
This project simulates how short DNA reads can be reconstructed into a full genome sequence using:
- k-mer generation
- Graph construction (De Bruijn–like graph)
- Eulerian path traversal
- Read simulation
- Error introduction
The graph operations are implemented using the R package igraph
## 🧠  **Background**
In de novo genome assembly:
1. A DNA sequence is broken into short fragments (reads).
2. Reads are decomposed into overlapping k-mers.
3. A graph is constructed where:
  -  Nodes represent k-mers
  -  Edges represent (k-1)-overlaps
4. If the graph satisfies Eulerian conditions, an Eulerian path reconstructs the original sequence.

This project simulates that full pipeline.

## **Functions**
### 1. kmerize(dna_sequences, k)
Generates all possible k-mers from a DNA sequence.

Parameters:
- dna_sequences – DNA string
- k – k-mer size

Returns:

Vector of k-mers

### 2. kmerGraph(kmers)
Constructs a graph where edges represent (k-1) overlaps between k-mers.

Returns:

An igraph graph object

### 3. isEulerian(graph)

If the graph is Eulerian, computes the Eulerian path and reconstructs the DNA sequence.

Returns:

Assembled DNA sequence or NULL

### 4. getAssembly(graph)

If the graph is Eulerian, computes the Eulerian path and reconstructs the DNA sequence.

Returns:

Assembled DNA sequence or NULL

### 5. readWriter(dna_sequence, readSize, numReads)

Simulates random sequencing reads from a DNA sequence.

Parameters:
- dna_sequence – original DNA
- readSize – length of each read
- numReads – number of reads to generate

Returns:

Vector of reads

### 6. readWronger(reads, errorRate)

Introduces random base substitution errors into reads.

Parameters:
- reads – vector of reads
- errorRate – proportion of bases to mutate (0–1)

Returns:

Modified reads with simulated sequencing errors

# What This Project Demonstrates
- Core principles of graph-based genome assembly
- How Eulerian paths enable sequence reconstruction
- Impact of read sampling and sequencing errors
- Practical use of graph theory in bioinformatics

# Limitations
- Does not implement a full De Bruijn graph with prefix/suffix nodes
- No coverage modeling
- No error correction algorithm
- Not optimized for large genomes
- Eulerian check assumes undirected even-degree condition
  
This is an educational simulation rather than a production-level assembler.

# Future Improvements: 
- Implement true De Bruijn graph (k-1 nodes)
- Add visualization of graph structure
- Add error correction module
- Handle non-Eulerian graphs
- Optimize performance for larger sequences
