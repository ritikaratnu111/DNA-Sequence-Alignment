# DNA-Sequence-Alignment

Genome assembly is a key problem in the area of computational genomics.Our assembly algorithm has been designed with the following key objectives:
1. The algorithm must be parallel and scalable for easy deployment on GPUs and FPGAs for acceleration.
2. The algorithm must preserve mutations of the genome.
We use sequence alignment to group reads into clusters. Clustering incurs one-time pre-processing cost. De-novo assembly is performed on a reduced read set, and a lot of unnecessary comparisons are prevented. Reads that align to locations within a distance or radius r form a cluster.  Value of r can be decided based on the desired parallelism and memory resource constraints. 

ART simulator was used to generate synthetic reads. The simuator mimicks the error profiles of the popular sequencing platforms to generate sequencing reads.
ACANA (a utility provided in the ART package) aligner was used to generate SAM file. fasta_converter.cpp and sam_filter.cpp are the files used for this purpose. assesmN files are the contigs reported from ART simulator wher N is the number of cluster.

A hash-map algorithm is used for de-novo assembly in the files assembler.cpp and block_assembler.cpp.



