Here you will find all the code used in the analyses of Gonçalves et al. 2019 and the tutorials presented in Gonçalves et al. Available at [https://www.sciencedirect.com/science/article/pii/S1055790318306602] 


### I. genome_assembly

In this first pipeline, we present the code used to processes Illumina reads using bbtools software to remove adaptors, phix, quality trim and to assemble reads into contigs.

### II. phylogenetic_analysis

Here we aligned genes from the plastid genome using two software, MAFFT and MACSE.  Then the alignments were used as input in a series of phylogenetic analyses using Maximum Likelihood and Multispecies Coalescent methods.

### III. treespace  
Now we used two algorithms (Robinson-Foulds and Kendal-Colijn) to explore the tree landscape of species and gene tree topologies inferred in the previous step. We used Principal Coordinate Analyses to assess the dissimilarity between the tree topologies inferred with different methods.

### IV. phylogenetic_signal  
We detected that the phylogenies inferred with different methods were incongruent in five nodes at different taxonomic levels. For three alternative topologies, we used IQ-TREE to calculate the maximum likelihood support of each site and gene. We used this information to calculate the delta site-wise and gene-wise log-likelihood support of alternative topologies.  

References:  

__Data article:__  
Gonçalves, DJP, Simpson, BB, Shimizu, GH, Jansen, RK, Ortiz, EM. 2019. Genome assembly and phylogenomic data analyses using plastid data: contrasting species tree estimation methods. Data in Brief. 25, 1-6
https://doi.org/10.1016/j.dib.2019.104271 
  
__Research article:__  
Gonçalves, DJP, Shimizu, GH, Ortiz, EMV, Simpson, BB, Jansen, RK. 2019. Incongruence between species tree and gene trees and phylogenetic signal variation in plastid genes. Molecular Phylogenetics and Evolution. 138, 219-232 https://doi.org/10.1016/j.ympev.2019.05.022 
