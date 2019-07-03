***
## II. PHYLOGENETIC ANALYSIS
**Authors:** Deise JP Goncalves et al.  
**Date:** 2018/06/04  
**Citation:**   
### Goncalves, DJP., Simpson, BB, Ortiz, EM, Shimizu, GH, Jansen, RK. Incongruence between species tree and gene trees and phylogenetic signal variation in plastid genes. Molecular Phylogenetics and Evolution. In review  
### DJP Goncalves, BB Simpson, GH Shimizu, RK Jansen, EM Ortiz. Genome assembly and phylogenomic data analyses using plastid data: contrasting species tree estimation methods. Data in Brief. Submitted  
***

We analysed a dataset of 78 plastid genes of 94 species (91 rosids and 3 asterids as outgroup). We generated four datasets:  
* Genes (CDS + introns, nucleotides) - **ge**
* Exons (CDS, nucleotides) - **ex**
* Codon-aligned (CDS, nucleotides) - **cd**
* Amino Acid (CDS, translated) - **aa**

Genes, Exons, and Amino Acid datasets were aligned in MAFFT; Codon dataset was aligned in MACSE. All datasets were inspected and concatenated in Geneious.

***

### `MAFFT and MACSE`
Commands below were used for all 78 protein-coding genes from the plastome. Here is an example of the code used for the first three genes:
```
#MAFFT
mafft --thread 12 --maxiterate 1000 --op 1.0 --genafpair accD.fasta > aligned_accD.fasta
mafft --thread 12 --maxiterate 1000 --op 1.0 --genafpair atpA.fasta > aligned_atpA.fasta
mafft --thread 12 --maxiterate 1000 --op 1.0 --genafpair atpB.fasta > aligned_atpB.fasta

#MACSE
java -jar macse_v0.9b1.jar -i accD_CDS_.fasta -o ../Macse_Output/macse_aligned_accD_CDS_.fasta -d 11 -f -30 -g -7 -s -100 -x -1
java -jar macse_v0.9b1.jar -i atpA_CDS.fasta -o ../Macse_Output/macse_aligned_atpA_CDS.fasta -d 11 -f -30 -g -7 -s -100 -x -1
java -jar macse_v0.9b1.jar -i atpB_CDS.fasta -o ../Macse_Output/macse_aligned_atpB_CDS.fasta -d 11 -f -30 -g -7 -s -100 -x -1
```
Phylogenetic analysis were performed using a supercomputer (Texas Adavance Computing Center-TACC, Lonestar5). Here we present the code used in each step of the maximum likelihood inference in IQ-TREE and the inferences using multispecies coalescent method performed in SVDquartets and ASTRAL-II.  
After performing the alignments, phylip files were visually inspected and exported as single gene and a concatenated alignment using Geneious v. 7.1.9 (Kearse et al., 2012).

***

### `IQ-TREE`


Tree search using concatenated alignments:
```
#gene alignment
iqtree-omp -nt 12 -st DNA -s Concatenated_Genes.phy -m TEST -AICc -b 1000 -pre geconc_mf

#exon alignment
iqtree-omp -nt 12 -st DNA -s Concatenated_Exons.phy -m TEST -AICc -b 1000 -pre exconc_mf

#codon alignment
iqtree-omp -nt 12 -st DNA -s codon_concatenated.phy -m TEST -AICc -b 1000 -pre cdconc_mf

#amino acid alignment
iqtree-omp -nt 12 -st AA -s Concatenated_AA.phy -msub chloroplast -m TEST -AICc -b 1000 -pre aaconc_mf
```

ModelFinder and tree search commands used for the concatenated alignments with partition by gene:
```
#Model Finder for all 4 matrices
iqtree-omp -nt AUTO -st DNA -s Concatenated_Genes.phy -spp Genes_conc_part.txt -m TEST -AICc -pre gegepart_mf_only
iqtree-omp -nt AUTO -st DNA -s Concatenated_Exons.phy -spp Exons_conc_part.txt -m TEST -AICc -pre excgepart_mf_only
iqtree-omp -nt AUTO -st DNA -s codon_concatenated.phy -spp codon_conc_part.txt -m TEST -AICc -pre cdcgepart_mf_only
iqtree-omp -nt AUTO -st AA  -s Concatenated_AA.phy    -spp Aa_conc_part.txt    -m TEST -AICc -pre aagepart_mf_only

#The previous code generates a NEXUS file (e.g. cdcgepart_mf_only.best_scheme.nex) containing the best-fit models of substitution

#gene alignment
iqtree-omp -nt 9 -st DNA -s Concatenated_Genes.phy -spp geconc_mf_only.best_scheme.nex -b 1000 -pre gegepart_ts

#exon alignment
iqtree-omp -nt 9 -st DNA -s Concatenated_Exons.phy -spp exconc_mf_only.best_scheme.nex -b 1000 -pre exgepart_ts

#codon alignment
iqtree-omp -nt 9 -st DNA -s codon_concatenated.phy -spp cdconc_mf_only.best_scheme.nex -b 1000 -pre cdgepart_ts

#amino acid alignment
iqtree-omp -nt 9 -st AA -s Concatenated_AA.phy -spp aaconc_mf_only.best_scheme.nex -b 1000 -pre aagepart_ts
```

ModelFinder and tree search of concatenated alignment with best scheme partition:
```
#gene alignment
#model finder
iqtree-omp -nt AUTO -st DNA -s concatenated_ge.phy -spp ge_part.txt -m TESTMERGEONLY -AICc -pre ge_mf
#tree search
sleep 0 && iqtree-omp -nt 12 -st DNA -s concatenated_ge.phy -spp ge_mf.best_scheme.nex -b 250 -pre gemfts_search1
sleep 1 && iqtree-omp -nt 12 -st DNA -s concatenated_ge.phy -spp ge_mf.best_scheme.nex -b 250 -pre gemfts_search2
sleep 2 && iqtree-omp -nt 12 -st DNA -s concatenated_ge.phy -spp ge_mf.best_scheme.nex -b 250 -pre gemfts_search3
sleep 3 && iqtree-omp -nt 12 -st DNA -s concatenated_ge.phy -spp ge_mf.best_scheme.nex -b 250 -pre gemfts_search4

#exon alignment
#model finder
iqtree-omp -nt AUTO -st DNA -s concatenated_ex.phy -spp ex_part.txt -m TESTMERGEONLY -AICc -pre ex_mf
#tree search
sleep 0 && iqtree-omp -nt 12 -st DNA -s concatenated_ex.phy -spp ex_mf.best_scheme.nex -b 250 -pre exmfts_search1
sleep 1 && iqtree-omp -nt 12 -st DNA -s concatenated_ex.phy -spp ex_mf.best_scheme.nex -b 250 -pre exmfts_search2
sleep 2 && iqtree-omp -nt 12 -st DNA -s concatenated_ex.phy -spp ex_mf.best_scheme.nex -b 250 -pre exmfts_search3
sleep 3 && iqtree-omp -nt 12 -st DNA -s concatenated_ex.phy -spp ex_mf.best_scheme.nex -b 250 -pre exmfts_search4

#codon alignment *
#model finder
iqtree-omp -nt AUTO -st DNA -s cd_concatenated.phy -spp cd_part.txt -msub chloroplast -m TESTONLY -AICc -pre cd_mf
#tree search
sleep 0 && iqtree-omp -nt 16 -st DNA -s cd_concatenated.phy -spp cd_mf.best_scheme.nex -b 250 -pre cdmfts_search1
sleep 1 && iqtree-omp -nt 16 -st DNA -s cd_concatenated.phy -spp cd_mf.best_scheme.nex -b 250 -pre cdmfts_search2
sleep 2 && iqtree-omp -nt 16 -st DNA -s cd_concatenated.phy -spp cd_mf.best_scheme.nex -b 250 -pre cdmfts_search3
sleep 3 && iqtree-omp -nt 16 -st DNA -s cd_concatenated.phy -spp cd_mf.best_scheme.nex -b 250 -pre cdmfts_search4

#amino acid alignment
#model finder
iqtree-omp -nt AUTO -st AA -s concatenated_aa.phy -spp aa_part.txt -m TESTMERGEONLY -msub chloroplast -AICc -pre aa_mf
#tree search
sleep 0 && iqtree-omp -nt 12 -st AA -s concatenated_aa.phy -spp aa_mf.best_scheme.nex -b 250 -pre aamfts_search1
sleep 1 && iqtree-omp -nt 12 -st AA -s concatenated_aa.phy -spp aa_mf.best_scheme.nex -b 250 -pre aamfts_search2
sleep 2 && iqtree-omp -nt 12 -st AA -s concatenated_aa.phy -spp aa_mf.best_scheme.nex -b 250 -pre aamfts_search3
sleep 3 && iqtree-omp -nt 12 -st AA -s concatenated_aa.phy -spp aa_mf.best_scheme.nex -b 250 -pre aamfts_search4
```
\* For the codon alignment in ModelFinder step of gene partition analysis, the `cd_part.txt` file used as input file had genes coordinates.  

Tree search of individual gene alignments (using the same models inferred in best scheme partition analyses):
```
#gene alignment
iqtree-omp -nt 3 -st DNA -s accD_ge.phy  -m GTR+I+G4 -b 100 -pre accD_ge_ts
iqtree-omp -nt 3 -st DNA -s atpA_ge.phy  -m GTR+I+G4 -b 100 -pre atpA_ge_ts
iqtree-omp -nt 3 -st DNA -s atpB_ge.phy  -m GTR+I+G4 -b 100 -pre atpB_ge_ts
iqtree-omp -nt 3 -st DNA -s atpE_ge.phy  -m TVM+I+G4 -b 100 -pre atpE_ge_ts
iqtree-omp -nt 3 -st DNA -s atpF_ge.phy  -m GTR+I+G4 -b 100 -pre atpF_ge_ts
iqtree-omp -nt 3 -st DNA -s atpH_ge.phy  -m K2P+G4   -b 100 -pre atpH_ge_ts
iqtree-omp -nt 3 -st DNA -s atpI_ge.phy  -m TVM+I+G4 -b 100 -pre atpI_ge_ts
iqtree-omp -nt 3 -st DNA -s ccsA_ge.phy  -m TVM+I+G4 -b 100 -pre ccsA_ge_ts
iqtree-omp -nt 3 -st DNA -s cemA_ge.phy  -m TVM+I+G4 -b 100 -pre cemA_ge_ts
iqtree-omp -nt 3 -st DNA -s clpP_ge.phy  -m GTR+I+G4 -b 100 -pre clpP_ge_ts
iqtree-omp -nt 3 -st DNA -s matK_ge.phy  -m TVM+I+G4 -b 100 -pre matK_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhA_ge.phy  -m TVM+I+G4 -b 100 -pre ndhA_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhB_ge.phy  -m TVM+I+G4 -b 100 -pre ndhB_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhC_ge.phy  -m GTR+I+G4 -b 100 -pre ndhC_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhD_ge.phy  -m GTR+I+G4 -b 100 -pre ndhD_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhE_ge.phy  -m TVM+G4   -b 100 -pre ndhE_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhF_ge.phy  -m TVM+I+G4 -b 100 -pre ndhF_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhG_ge.phy  -m GTR+I+G4 -b 100 -pre ndhG_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhH_ge.phy  -m GTR+I+G4 -b 100 -pre ndhH_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhI_ge.phy  -m GTR+I+G4 -b 100 -pre ndhI_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhJ_ge.phy  -m TVM+I+G4 -b 100 -pre ndhJ_ge_ts
iqtree-omp -nt 3 -st DNA -s ndhK_ge.phy  -m GTR+I+G4 -b 100 -pre ndhK_ge_ts
iqtree-omp -nt 3 -st DNA -s petA_ge.phy  -m TVM+I+G4 -b 100 -pre petA_ge_ts
iqtree-omp -nt 3 -st DNA -s petB_ge.phy  -m TVM+I+G4 -b 100 -pre petB_ge_ts
iqtree-omp -nt 3 -st DNA -s petD_ge.phy  -m GTR+I+G4 -b 100 -pre petD_ge_ts
iqtree-omp -nt 3 -st DNA -s petG_ge.phy  -m JC       -b 100 -pre petG_ge_ts
iqtree-omp -nt 3 -st DNA -s petL_ge.phy  -m JC       -b 100 -pre petL_ge_ts
iqtree-omp -nt 3 -st DNA -s petN_ge.phy  -m JC       -b 100 -pre petN_ge_ts
iqtree-omp -nt 3 -st DNA -s psaA_ge.phy  -m GTR+I+G4 -b 100 -pre psaA_ge_ts
iqtree-omp -nt 3 -st DNA -s psaB_ge.phy  -m TVM+I+G4 -b 100 -pre psaB_ge_ts
iqtree-omp -nt 3 -st DNA -s psaC_ge.phy  -m K2P+I+G4 -b 100 -pre psaC_ge_ts
iqtree-omp -nt 3 -st DNA -s psaI_ge.phy  -m JC       -b 100 -pre psaI_ge_ts
iqtree-omp -nt 3 -st DNA -s psaJ_ge.phy  -m JC       -b 100 -pre psaJ_ge_ts
iqtree-omp -nt 3 -st DNA -s psbA_ge.phy  -m GTR+I+G4 -b 100 -pre psbA_ge_ts
iqtree-omp -nt 3 -st DNA -s psbB_ge.phy  -m TVM+I+G4 -b 100 -pre psbB_ge_ts
iqtree-omp -nt 3 -st DNA -s psbC_ge.phy  -m TVM+I+G4 -b 100 -pre psbC_ge_ts
iqtree-omp -nt 3 -st DNA -s psbD_ge.phy  -m TVM+I+G4 -b 100 -pre psbD_ge_ts
iqtree-omp -nt 3 -st DNA -s psbE_ge.phy  -m K2P+I+G4 -b 100 -pre psbE_ge_ts
iqtree-omp -nt 3 -st DNA -s psbF_ge.phy  -m JC       -b 100 -pre psbF_ge_ts
iqtree-omp -nt 3 -st DNA -s psbH_ge.phy  -m K2P+G4   -b 100 -pre psbH_ge_ts
iqtree-omp -nt 3 -st DNA -s psbI_ge.phy  -m JC       -b 100 -pre psbI_ge_ts
iqtree-omp -nt 3 -st DNA -s psbJ_ge.phy  -m JC       -b 100 -pre psbJ_ge_ts
iqtree-omp -nt 3 -st DNA -s psbK_ge.phy  -m JC+G4    -b 100 -pre psbK_ge_ts
iqtree-omp -nt 3 -st DNA -s psbL_ge.phy  -m JC       -b 100 -pre psbL_ge_ts
iqtree-omp -nt 3 -st DNA -s psbM_ge.phy  -m JC       -b 100 -pre psbM_ge_ts
iqtree-omp -nt 3 -st DNA -s psbN_ge.phy  -m JC       -b 100 -pre psbN_ge_ts
iqtree-omp -nt 3 -st DNA -s psbT_ge.phy  -m JC       -b 100 -pre psbT_ge_ts
iqtree-omp -nt 3 -st DNA -s psbZ_ge.phy  -m JC       -b 100 -pre psbZ_ge_ts
iqtree-omp -nt 3 -st DNA -s rbcL_ge.phy  -m GTR+I+G4 -b 100 -pre rbcL_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl14_ge.phy -m K3Pu+G4  -b 100 -pre rpl14_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl16_ge.phy -m TVM+I+G4 -b 100 -pre rpl16_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl2_ge.phy  -m TVM+I+G4 -b 100 -pre rpl2_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl20_ge.phy -m TVM+G4   -b 100 -pre rpl20_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl22_ge.phy -m TVM+I+G4 -b 100 -pre rpl22_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl23_ge.phy -m TVM+G4   -b 100 -pre rpl23_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl32_ge.phy -m TVM+G4   -b 100 -pre rpl32_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl33_ge.phy -m TPM3u+G4 -b 100 -pre rpl33_ge_ts
iqtree-omp -nt 3 -st DNA -s rpl36_ge.phy -m JC       -b 100 -pre rpl36_ge_ts
iqtree-omp -nt 3 -st DNA -s rpoA_ge.phy  -m TVM+I+G4 -b 100 -pre rpoA_ge_ts
iqtree-omp -nt 3 -st DNA -s rpoB_ge.phy  -m GTR+I+G4 -b 100 -pre rpoB_ge_ts
iqtree-omp -nt 3 -st DNA -s rpoC1_ge.phy -m TVM+I+G4 -b 100 -pre rpoC1_ge_ts
iqtree-omp -nt 3 -st DNA -s rpoC2_ge.phy -m GTR+I+G4 -b 100 -pre rpoC2_ge_ts
iqtree-omp -nt 3 -st DNA -s rps11_ge.phy -m TVM+G4   -b 100 -pre rps11_ge_ts
iqtree-omp -nt 3 -st DNA -s rps12_ge.phy -m GTR+I+G4 -b 100 -pre rps12_ge_ts
iqtree-omp -nt 3 -st DNA -s rps14_ge.phy -m TVM+G4   -b 100 -pre rps14_ge_ts
iqtree-omp -nt 3 -st DNA -s rps15_ge.phy -m K3Pu+G4  -b 100 -pre rps15_ge_ts
iqtree-omp -nt 3 -st DNA -s rps16_ge.phy -m GTR+I+G4 -b 100 -pre rps16_ge_ts
iqtree-omp -nt 3 -st DNA -s rps18_ge.phy -m TVM+I+G4 -b 100 -pre rps18_ge_ts
iqtree-omp -nt 3 -st DNA -s rps19_ge.phy -m K3Pu+G4  -b 100 -pre rps19_ge_ts
iqtree-omp -nt 3 -st DNA -s rps2_ge.phy  -m GTR+I+G4 -b 100 -pre rps2_ge_ts
iqtree-omp -nt 3 -st DNA -s rps3_ge.phy  -m TVM+I+G4 -b 100 -pre rps3_ge_ts
iqtree-omp -nt 3 -st DNA -s rps4_ge.phy  -m TVM+G4   -b 100 -pre rps4_ge_ts
iqtree-omp -nt 3 -st DNA -s rps7_ge.phy  -m TVM+G4   -b 100 -pre rps7_ge_ts
iqtree-omp -nt 3 -st DNA -s rps8_ge.phy  -m TVM+G4   -b 100 -pre rps8_ge_ts
iqtree-omp -nt 3 -st DNA -s ycf1_ge.phy  -m GTR+I+G4 -b 100 -pre ycf1_ge_ts
iqtree-omp -nt 3 -st DNA -s ycf2_ge.phy  -m GTR+G4   -b 100 -pre ycf2_ge_ts
iqtree-omp -nt 3 -st DNA -s ycf3_ge.phy  -m GTR+I+G4 -b 100 -pre ycf3_ge_ts
iqtree-omp -nt 3 -st DNA -s ycf4_ge.phy  -m TVM+G4   -b 100 -pre ycf4_ge_ts

#exon alignment
iqtree-omp -nt 3 -st DNA -s accD_ex.phy  -m GTR+R5   -b 100 -pre accD_ex_ts
iqtree-omp -nt 3 -st DNA -s atpA_ex.phy  -m GTR+R5   -b 100 -pre atpA_ex_ts
iqtree-omp -nt 3 -st DNA -s atpB_ex.phy  -m GTR+R4   -b 100 -pre atpB_ex_ts
iqtree-omp -nt 3 -st DNA -s atpE_ex.phy  -m TVM+R3   -b 100 -pre atpE_ex_ts
iqtree-omp -nt 3 -st DNA -s atpF_ex.phy  -m TVM+R4   -b 100 -pre atpF_ex_ts
iqtree-omp -nt 3 -st DNA -s atpH_ex.phy  -m K2P+G4   -b 100 -pre atpH_ex_ts
iqtree-omp -nt 3 -st DNA -s atpI_ex.phy  -m TVM+I+G4 -b 100 -pre atpI_ex_ts
iqtree-omp -nt 3 -st DNA -s ccsA_ex.phy  -m GTR+R5   -b 100 -pre ccsA_ex_ts
iqtree-omp -nt 3 -st DNA -s cemA_ex.phy  -m TVM+R5   -b 100 -pre cemA_ex_ts
iqtree-omp -nt 3 -st DNA -s clpP_ex.phy  -m GTR+R4   -b 100 -pre clpP_ex_ts
iqtree-omp -nt 3 -st DNA -s matK_ex.phy  -m TVM+R5   -b 100 -pre matK_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhA_ex.phy  -m TVM+R4   -b 100 -pre ndhA_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhB_ex.phy  -m GTR+R3   -b 100 -pre ndhB_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhC_ex.phy  -m TVM+R3   -b 100 -pre ndhC_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhD_ex.phy  -m TVM+R5   -b 100 -pre ndhD_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhE_ex.phy  -m TVM+G4   -b 100 -pre ndhE_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhF_ex.phy  -m TVM+R6   -b 100 -pre ndhF_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhG_ex.phy  -m GTR+I+G4 -b 100 -pre ndhG_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhH_ex.phy  -m TVM+R4   -b 100 -pre ndhH_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhI_ex.phy  -m GTR+R4   -b 100 -pre ndhI_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhJ_ex.phy  -m TVM+R3   -b 100 -pre ndhJ_ex_ts
iqtree-omp -nt 3 -st DNA -s ndhK_ex.phy  -m TVM+R4   -b 100 -pre ndhK_ex_ts
iqtree-omp -nt 3 -st DNA -s petA_ex.phy  -m TVM+R4   -b 100 -pre petA_ex_ts
iqtree-omp -nt 3 -st DNA -s petB_ex.phy  -m TVM+R4   -b 100 -pre petB_ex_ts
iqtree-omp -nt 3 -st DNA -s petD_ex.phy  -m GTR+I+G4 -b 100 -pre petD_ex_ts
iqtree-omp -nt 3 -st DNA -s petG_ex.phy  -m JC       -b 100 -pre petG_ex_ts
iqtree-omp -nt 3 -st DNA -s petL_ex.phy  -m JC       -b 100 -pre petL_ex_ts
iqtree-omp -nt 3 -st DNA -s petN_ex.phy  -m JC       -b 100 -pre petN_ex_ts
iqtree-omp -nt 3 -st DNA -s psaA_ex.phy  -m GTR+R4   -b 100 -pre psaA_ex_ts
iqtree-omp -nt 3 -st DNA -s psaB_ex.phy  -m TVM+R4   -b 100 -pre psaB_ex_ts
iqtree-omp -nt 3 -st DNA -s psaC_ex.phy  -m K2P+I+G4 -b 100 -pre psaC_ex_ts
iqtree-omp -nt 3 -st DNA -s psaI_ex.phy  -m JC       -b 100 -pre psaI_ex_ts
iqtree-omp -nt 3 -st DNA -s psaJ_ex.phy  -m JC       -b 100 -pre psaJ_ex_ts
iqtree-omp -nt 3 -st DNA -s psbA_ex.phy  -m GTR+R4   -b 100 -pre psbA_ex_ts
iqtree-omp -nt 3 -st DNA -s psbB_ex.phy  -m GTR+R4   -b 100 -pre psbB_ex_ts
iqtree-omp -nt 3 -st DNA -s psbC_ex.phy  -m GTR+R4   -b 100 -pre psbC_ex_ts
iqtree-omp -nt 3 -st DNA -s psbD_ex.phy  -m TVM+R4   -b 100 -pre psbD_ex_ts
iqtree-omp -nt 3 -st DNA -s psbE_ex.phy  -m K2P+I+G4 -b 100 -pre psbE_ex_ts
iqtree-omp -nt 3 -st DNA -s psbF_ex.phy  -m JC       -b 100 -pre psbF_ex_ts
iqtree-omp -nt 3 -st DNA -s psbH_ex.phy  -m K2P+G4   -b 100 -pre psbH_ex_ts
iqtree-omp -nt 3 -st DNA -s psbI_ex.phy  -m JC       -b 100 -pre psbI_ex_ts
iqtree-omp -nt 3 -st DNA -s psbJ_ex.phy  -m JC       -b 100 -pre psbJ_ex_ts
iqtree-omp -nt 3 -st DNA -s psbK_ex.phy  -m JC+G4    -b 100 -pre psbK_ex_ts
iqtree-omp -nt 3 -st DNA -s psbL_ex.phy  -m JC       -b 100 -pre psbL_ex_ts
iqtree-omp -nt 3 -st DNA -s psbM_ex.phy  -m JC       -b 100 -pre psbM_ex_ts
iqtree-omp -nt 3 -st DNA -s psbN_ex.phy  -m JC       -b 100 -pre psbN_ex_ts
iqtree-omp -nt 3 -st DNA -s psbT_ex.phy  -m JC       -b 100 -pre psbT_ex_ts
iqtree-omp -nt 3 -st DNA -s psbZ_ex.phy  -m JC       -b 100 -pre psbZ_ex_ts
iqtree-omp -nt 3 -st DNA -s rbcL_ex.phy  -m GTR+R5   -b 100 -pre rbcL_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl14_ex.phy -m K3Pu+G4  -b 100 -pre rpl14_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl16_ex.phy -m K3Pu+R3  -b 100 -pre rpl16_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl20_ex.phy -m TVM+R4   -b 100 -pre rpl20_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl22_ex.phy -m TVM+R5   -b 100 -pre rpl22_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl23_ex.phy -m TVM+R5   -b 100 -pre rpl23_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl2_ex.phy  -m GTR+R5   -b 100 -pre rpl2_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl32_ex.phy -m TVM+G4   -b 100 -pre rpl32_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl33_ex.phy -m TPM3+G4  -b 100 -pre rpl33_ex_ts
iqtree-omp -nt 3 -st DNA -s rpl36_ex.phy -m JC       -b 100 -pre rpl36_ex_ts
iqtree-omp -nt 3 -st DNA -s rpoA_ex.phy  -m TVM+R4   -b 100 -pre rpoA_ex_ts
iqtree-omp -nt 3 -st DNA -s rpoB_ex.phy  -m GTR+R5   -b 100 -pre rpoB_ex_ts
iqtree-omp -nt 3 -st DNA -s rpoC1_ex.phy -m TVM+R5   -b 100 -pre rpoC1_ex_ts
iqtree-omp -nt 3 -st DNA -s rpoC2_ex.phy -m GTR+R5   -b 100 -pre rpoC2_ex_ts
iqtree-omp -nt 3 -st DNA -s rps11_ex.phy -m TVM+R3   -b 100 -pre rps11_ex_ts
iqtree-omp -nt 3 -st DNA -s rps12_ex.phy -m TVM+G4   -b 100 -pre rps12_ex_ts
iqtree-omp -nt 3 -st DNA -s rps14_ex.phy -m TVM+G4   -b 100 -pre rps14_ex_ts
iqtree-omp -nt 3 -st DNA -s rps15_ex.phy -m K3Pu+G4  -b 100 -pre rps15_ex_ts
iqtree-omp -nt 3 -st DNA -s rps16_ex.phy -m TVM+G4   -b 100 -pre rps16_ex_ts
iqtree-omp -nt 3 -st DNA -s rps18_ex.phy -m TVM+R4   -b 100 -pre rps18_ex_ts
iqtree-omp -nt 3 -st DNA -s rps19_ex.phy -m TVM+R2   -b 100 -pre rps19_ex_ts
iqtree-omp -nt 3 -st DNA -s rps2_ex.phy  -m GTR+R4   -b 100 -pre rps2_ex_ts
iqtree-omp -nt 3 -st DNA -s rps3_ex.phy  -m TVM+R4   -b 100 -pre rps3_ex_ts
iqtree-omp -nt 3 -st DNA -s rps4_ex.phy  -m GTR+R3   -b 100 -pre rps4_ex_ts
iqtree-omp -nt 3 -st DNA -s rps7_ex.phy  -m TVM+R3   -b 100 -pre rps7_ex_ts
iqtree-omp -nt 3 -st DNA -s rps8_ex.phy  -m TVM+G4   -b 100 -pre rps8_ex_ts
iqtree-omp -nt 3 -st DNA -s ycf1_ex.phy  -m GTR+R5   -b 100 -pre ycf1_ex_ts
iqtree-omp -nt 3 -st DNA -s ycf2_ex.phy  -m GTR+R4   -b 100 -pre ycf2_ex_ts
iqtree-omp -nt 3 -st DNA -s ycf3_ex.phy  -m TVM+R3   -b 100 -pre ycf3_ex_ts
iqtree-omp -nt 3 -st DNA -s ycf4_ex.phy  -m TVM+R4   -b 100 -pre ycf4_ex_ts

#codon alignment
iqtree-omp -nt 3 -st DNA -s accD_cd.phy  -m GTR+I+G4  -b 100 -pre accD_cdts
iqtree-omp -nt 3 -st DNA -s atpA_cd.phy  -m GTR+I+G4  -b 100 -pre atpA_cdts
iqtree-omp -nt 3 -st DNA -s atpB_cd.phy  -m GTR+I+G4  -b 100 -pre atpB_cdts
iqtree-omp -nt 3 -st DNA -s atpE_cd.phy  -m TVM+I+G4  -b 100 -pre atpE_cdts
iqtree-omp -nt 3 -st DNA -s atpF_cd.phy  -m TVM+I+G4  -b 100 -pre atpF_cdts
iqtree-omp -nt 3 -st DNA -s atpH_cd.phy  -m K2P+G4    -b 100 -pre atpH_cdts
iqtree-omp -nt 3 -st DNA -s atpI_cd.phy  -m TVM+I+G4  -b 100 -pre atpI_cdts
iqtree-omp -nt 3 -st DNA -s ccsA_cd.phy  -m TVM+I+G4  -b 100 -pre ccsA_cdts
iqtree-omp -nt 3 -st DNA -s cemA_cd.phy  -m TVM+I+G4  -b 100 -pre cemA_cdts
iqtree-omp -nt 3 -st DNA -s clpP_cd.phy  -m TVM+G4    -b 100 -pre clpP_cdts
iqtree-omp -nt 3 -st DNA -s matK_cd.phy  -m TVM+I+G4  -b 100 -pre matK_cdts
iqtree-omp -nt 3 -st DNA -s ndhA_cd.phy  -m GTR+I+G4  -b 100 -pre ndhA_cdts
iqtree-omp -nt 3 -st DNA -s ndhB_cd.phy  -m TVM+I+G4  -b 100 -pre ndhB_cdts
iqtree-omp -nt 3 -st DNA -s ndhC_cd.phy  -m GTR+I+G4  -b 100 -pre ndhC_cdts
iqtree-omp -nt 3 -st DNA -s ndhD_cd.phy  -m GTR+I+G4  -b 100 -pre ndhD_cdts
iqtree-omp -nt 3 -st DNA -s ndhE_cd.phy  -m TVM+G4    -b 100 -pre ndhE_cdts
iqtree-omp -nt 3 -st DNA -s ndhF_cd.phy  -m TVM+I+G4  -b 100 -pre ndhF_cdts
iqtree-omp -nt 3 -st DNA -s ndhG_cd.phy  -m GTR+I+G4  -b 100 -pre ndhG_cdts
iqtree-omp -nt 3 -st DNA -s ndhH_cd.phy  -m GTR+I+G4  -b 100 -pre ndhH_cdts
iqtree-omp -nt 3 -st DNA -s ndhI_cd.phy  -m TVM+I+G4  -b 100 -pre ndhI_cdts
iqtree-omp -nt 3 -st DNA -s ndhJ_cd.phy  -m TVM+I+G4  -b 100 -pre ndhJ_cdts
iqtree-omp -nt 3 -st DNA -s ndhK_cd.phy  -m GTR+I+G4  -b 100 -pre ndhK_cdts
iqtree-omp -nt 3 -st DNA -s petA_cd.phy  -m TVM+I+G4  -b 100 -pre petA_cdts
iqtree-omp -nt 3 -st DNA -s petB_cd.phy  -m TVM+I+G4  -b 100 -pre petB_cdts
iqtree-omp -nt 3 -st DNA -s petD_cd.phy  -m GTR+I+G4  -b 100 -pre petD_cdts
iqtree-omp -nt 3 -st DNA -s petG_cd.phy  -m JC        -b 100 -pre petG_cdts
iqtree-omp -nt 3 -st DNA -s petL_cd.phy  -m JC        -b 100 -pre petL_cdts
iqtree-omp -nt 3 -st DNA -s petN_cd.phy  -m JC        -b 100 -pre petN_cdts
iqtree-omp -nt 3 -st DNA -s psaA_cd.phy  -m GTR+I+G4  -b 100 -pre psaA_cdts
iqtree-omp -nt 3 -st DNA -s psaB_cd.phy  -m TVM+I+G4  -b 100 -pre psaB_cdts
iqtree-omp -nt 3 -st DNA -s psaC_cd.phy  -m K2P+I+G4  -b 100 -pre psaC_cdts
iqtree-omp -nt 3 -st DNA -s psaI_cd.phy  -m JC        -b 100 -pre psaI_cdts
iqtree-omp -nt 3 -st DNA -s psaJ_cd.phy  -m JC        -b 100 -pre psaJ_cdts
iqtree-omp -nt 3 -st DNA -s psbA_cd.phy  -m GTR+I+G4  -b 100 -pre psbA_cdts
iqtree-omp -nt 3 -st DNA -s psbB_cd.phy  -m TVM+I+G4  -b 100 -pre psbB_cdts
iqtree-omp -nt 3 -st DNA -s psbC_cd.phy  -m TVM+I+G4  -b 100 -pre psbC_cdts
iqtree-omp -nt 3 -st DNA -s psbD_cd.phy  -m TVM+I+G4  -b 100 -pre psbD_cdts
iqtree-omp -nt 3 -st DNA -s psbE_cd.phy  -m K2P+I+G4  -b 100 -pre psbE_cdts
iqtree-omp -nt 3 -st DNA -s psbF_cd.phy  -m JC        -b 100 -pre psbF_cdts
iqtree-omp -nt 3 -st DNA -s psbH_cd.phy  -m K2P+G4    -b 100 -pre psbH_cdts
iqtree-omp -nt 3 -st DNA -s psbI_cd.phy  -m JC        -b 100 -pre psbI_cdts
iqtree-omp -nt 3 -st DNA -s psbJ_cd.phy  -m JC        -b 100 -pre psbJ_cdts
iqtree-omp -nt 3 -st DNA -s psbK_cd.phy  -m JC        -b 100 -pre psbK_cdts
iqtree-omp -nt 3 -st DNA -s psbL_cd.phy  -m JC        -b 100 -pre psbL_cdts
iqtree-omp -nt 3 -st DNA -s psbM_cd.phy  -m JC        -b 100 -pre psbM_cdts
iqtree-omp -nt 3 -st DNA -s psbN_cd.phy  -m JC        -b 100 -pre psbN_cdts
iqtree-omp -nt 3 -st DNA -s psbT_cd.phy  -m JC        -b 100 -pre psbT_cdts
iqtree-omp -nt 3 -st DNA -s psbZ_cd.phy  -m JC        -b 100 -pre psbZ_cdts
iqtree-omp -nt 3 -st DNA -s rbcL_cd.phy  -m GTR+I+G4  -b 100 -pre rbcL_cdts
iqtree-omp -nt 3 -st DNA -s rpl14_cd.phy -m K3Pu+G4   -b 100 -pre rpl14_cdts
iqtree-omp -nt 3 -st DNA -s rpl16_cd.phy -m K3Pu+I+G4 -b 100 -pre rpl16_cdts
iqtree-omp -nt 3 -st DNA -s rpl20_cd.phy -m TVM+G4    -b 100 -pre rpl20_cdts
iqtree-omp -nt 3 -st DNA -s rpl22_cd.phy -m TVM+I+G4  -b 100 -pre rpl22_cdts
iqtree-omp -nt 3 -st DNA -s rpl23_cd.phy -m TVM+G4    -b 100 -pre rpl23_cdts
iqtree-omp -nt 3 -st DNA -s rpl2_cd.phy  -m TVM+I+G4  -b 100 -pre rpl2_cdts
iqtree-omp -nt 3 -st DNA -s rpl32_cd.phy -m TVM+G4    -b 100 -pre rpl32_cdts
iqtree-omp -nt 3 -st DNA -s rpl33_cd.phy -m K2P+G4    -b 100 -pre rpl33_cdts
iqtree-omp -nt 3 -st DNA -s rpl36_cd.phy -m JC        -b 100 -pre rpl36_cdts
iqtree-omp -nt 3 -st DNA -s rpoA_cd.phy  -m GTR+I+G4  -b 100 -pre rpoA_cdts
iqtree-omp -nt 3 -st DNA -s rpoB_cd.phy  -m GTR+I+G4  -b 100 -pre rpoB_cdts
iqtree-omp -nt 3 -st DNA -s rpoC1_cd.phy -m TVM+I+G4  -b 100 -pre rpoC1_cdts
iqtree-omp -nt 3 -st DNA -s rpoC2_cd.phy -m GTR+I+G4  -b 100 -pre rpoC2_cdts
iqtree-omp -nt 3 -st DNA -s rps11_cd.phy -m GTR+G4    -b 100 -pre rps11_cdts
iqtree-omp -nt 3 -st DNA -s rps12_cd.phy -m TVM+G4    -b 100 -pre rps12_cdts
iqtree-omp -nt 3 -st DNA -s rps14_cd.phy -m TVM+G4    -b 100 -pre rps14_cdts
iqtree-omp -nt 3 -st DNA -s rps15_cd.phy -m K3Pu+G4   -b 100 -pre rps15_cdts
iqtree-omp -nt 3 -st DNA -s rps16_cd.phy -m TVM+G4    -b 100 -pre rps16_cdts
iqtree-omp -nt 3 -st DNA -s rps18_cd.phy -m TVM+I+G4  -b 100 -pre rps18_cdts
iqtree-omp -nt 3 -st DNA -s rps19_cd.phy -m TVM+G4    -b 100 -pre rps19_cdts
iqtree-omp -nt 3 -st DNA -s rps2_cd.phy  -m TVM+I+G4  -b 100 -pre rps2_cdts
iqtree-omp -nt 3 -st DNA -s rps3_cd.phy  -m TVM+I+G4  -b 100 -pre rps3_cdts
iqtree-omp -nt 3 -st DNA -s rps4_cd.phy  -m TVM+G4    -b 100 -pre rps4_cdts
iqtree-omp -nt 3 -st DNA -s rps7_cd.phy  -m K3Pu+G4   -b 100 -pre rps7_cdts
iqtree-omp -nt 3 -st DNA -s rps8_cd.phy  -m TVM+G4    -b 100 -pre rps8_cdts
iqtree-omp -nt 3 -st DNA -s ycf1_cd.phy  -m GTR+I+G4  -b 100 -pre ycf1_cdts
iqtree-omp -nt 3 -st DNA -s ycf2_cd.phy  -m GTR+G4    -b 100 -pre ycf2_cdts
iqtree-omp -nt 3 -st DNA -s ycf3_cd.phy  -m TVM+I+G4  -b 100 -pre ycf3_cdts
iqtree-omp -nt 3 -st DNA -s ycf4_cd.phy  -m TVM+G4    -b 100 -pre ycf4_cdts

#amino acid alignment
iqtree-omp -nt 3 -st AA -s  accD_aa.phy  -m cpREV+F+I+G4 -b 100 -pre accD_aa_ts
iqtree-omp -nt 3 -st AA -s  atpA_aa.phy  -m cpREV+I+G4   -b 100 -pre atpA_aa_ts
iqtree-omp -nt 3 -st AA -s  atpB_aa.phy  -m cpREV+I+G4   -b 100 -pre atpB_aa_ts
iqtree-omp -nt 3 -st AA -s  atpE_aa.phy  -m cpREV        -b 100 -pre atpE_aa_ts
iqtree-omp -nt 3 -st AA -s  atpF_aa.phy  -m cpREV+G4     -b 100 -pre atpF_aa_ts
iqtree-omp -nt 3 -st AA -s  atpH_aa.phy  -m cpREV        -b 100 -pre atpH_aa_ts
iqtree-omp -nt 3 -st AA -s  atpI_aa.phy  -m cpREV+G4     -b 100 -pre atpI_aa_ts
iqtree-omp -nt 3 -st AA -s  ccsA_aa.phy  -m cpREV+I+G4   -b 100 -pre ccsA_aa_ts
iqtree-omp -nt 3 -st AA -s  cemA_aa.phy  -m cpREV+G4     -b 100 -pre cemA_aa_ts
iqtree-omp -nt 3 -st AA -s  clpP_aa.phy  -m cpREV+G4     -b 100 -pre clpP_aa_ts
iqtree-omp -nt 3 -st AA -s  matK_aa.phy  -m cpREV+F+I+G4 -b 100 -pre matK_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhA_aa.phy  -m cpREV+I+G4   -b 100 -pre ndhA_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhB_aa.phy  -m cpREV+F+I+G4 -b 100 -pre ndhB_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhC_aa.phy  -m cpREV        -b 100 -pre ndhC_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhD_aa.phy  -m cpREV+F+I+G4 -b 100 -pre ndhD_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhE_aa.phy  -m cpREV        -b 100 -pre ndhE_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhF_aa.phy  -m cpREV+F+I+G4 -b 100 -pre ndhF_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhG_aa.phy  -m cpREV        -b 100 -pre ndhG_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhH_aa.phy  -m cpREV+I+G4   -b 100 -pre ndhH_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhI_aa.phy  -m cpREV        -b 100 -pre ndhI_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhJ_aa.phy  -m cpREV        -b 100 -pre ndhJ_aa_ts
iqtree-omp -nt 3 -st AA -s  ndhK_aa.phy  -m cpREV+G4     -b 100 -pre ndhK_aa_ts
iqtree-omp -nt 3 -st AA -s  petA_aa.phy  -m cpREV+G4     -b 100 -pre petA_aa_ts
iqtree-omp -nt 3 -st AA -s  petB_aa.phy  -m cpREV+I+G4   -b 100 -pre petB_aa_ts
iqtree-omp -nt 3 -st AA -s  petD_aa.phy  -m cpREV        -b 100 -pre petD_aa_ts
iqtree-omp -nt 3 -st AA -s  petG_aa.phy  -m cpREV        -b 100 -pre petG_aa_ts
iqtree-omp -nt 3 -st AA -s  petL_aa.phy  -m cpREV        -b 100 -pre petL_aa_ts
iqtree-omp -nt 3 -st AA -s  petN_aa.phy  -m cpREV        -b 100 -pre petN_aa_ts
iqtree-omp -nt 3 -st AA -s  psaA_aa.phy  -m cpREV+I+G4   -b 100 -pre psaA_aa_ts
iqtree-omp -nt 3 -st AA -s  psaB_aa.phy  -m cpREV+F+I+G4 -b 100 -pre psaB_aa_ts
iqtree-omp -nt 3 -st AA -s  psaC_aa.phy  -m cpREV        -b 100 -pre psaC_aa_ts
iqtree-omp -nt 3 -st AA -s  psaI_aa.phy  -m cpREV        -b 100 -pre psaI_aa_ts
iqtree-omp -nt 3 -st AA -s  psaJ_aa.phy  -m cpREV        -b 100 -pre psaJ_aa_ts
iqtree-omp -nt 3 -st AA -s  psbA_aa.phy  -m cpREV+I+G4   -b 100 -pre psbA_aa_ts
iqtree-omp -nt 3 -st AA -s  psbB_aa.phy  -m cpREV+I+G4   -b 100 -pre psbB_aa_ts
iqtree-omp -nt 3 -st AA -s  psbC_aa.phy  -m cpREV+I+G4   -b 100 -pre psbC_aa_ts
iqtree-omp -nt 3 -st AA -s  psbD_aa.phy  -m cpREV+I+G4   -b 100 -pre psbD_aa_ts
iqtree-omp -nt 3 -st AA -s  psbE_aa.phy  -m cpREV        -b 100 -pre psbE_aa_ts
iqtree-omp -nt 3 -st AA -s  psbF_aa.phy  -m cpREV        -b 100 -pre psbF_aa_ts
iqtree-omp -nt 3 -st AA -s  psbH_aa.phy  -m cpREV        -b 100 -pre psbH_aa_ts
iqtree-omp -nt 3 -st AA -s  psbI_aa.phy  -m cpREV        -b 100 -pre psbI_aa_ts
iqtree-omp -nt 3 -st AA -s  psbJ_aa.phy  -m cpREV        -b 100 -pre psbJ_aa_ts
iqtree-omp -nt 3 -st AA -s  psbK_aa.phy  -m cpREV        -b 100 -pre psbK_aa_ts
iqtree-omp -nt 3 -st AA -s  psbL_aa.phy  -m cpREV        -b 100 -pre psbL_aa_ts
iqtree-omp -nt 3 -st AA -s  psbM_aa.phy  -m cpREV        -b 100 -pre psbM_aa_ts
iqtree-omp -nt 3 -st AA -s  psbN_aa.phy  -m cpREV        -b 100 -pre psbN_aa_ts
iqtree-omp -nt 3 -st AA -s  psbT_aa.phy  -m cpREV        -b 100 -pre psbT_aa_ts
iqtree-omp -nt 3 -st AA -s  psbZ_aa.phy  -m cpREV        -b 100 -pre psbZ_aa_ts
iqtree-omp -nt 3 -st AA -s  rbcL_aa.phy  -m cpREV+I+G4   -b 100 -pre rbcL_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl14_aa.phy -m cpREV        -b 100 -pre rpl14_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl16_aa.phy -m cpREV        -b 100 -pre rpl16_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl2_aa.phy  -m cpREV+G4     -b 100 -pre rpl2_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl20_aa.phy -m cpREV        -b 100 -pre rpl20_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl22_aa.phy -m cpREV+G4     -b 100 -pre rpl22_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl23_aa.phy -m cpREV        -b 100 -pre rpl23_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl32_aa.phy -m cpREV        -b 100 -pre rpl32_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl33_aa.phy -m cpREV        -b 100 -pre rpl33_aa_ts
iqtree-omp -nt 3 -st AA -s  rpl36_aa.phy -m cpREV        -b 100 -pre rpl36_aa_ts
iqtree-omp -nt 3 -st AA -s  rpoA_aa.phy  -m cpREV+G4     -b 100 -pre rpoA_aa_ts
iqtree-omp -nt 3 -st AA -s  rpoB_aa.phy  -m cpREV+I+G4   -b 100 -pre rpoB_aa_ts
iqtree-omp -nt 3 -st AA -s  rpoC1_aa.phy -m cpREV+G4     -b 100 -pre rpoC1_aa_ts
iqtree-omp -nt 3 -st AA -s  rpoC2_aa.phy -m cpREV+I+G4   -b 100 -pre rpoC2_aa_ts
iqtree-omp -nt 3 -st AA -s  rps11_aa.phy -m cpREV        -b 100 -pre rps11_aa_ts
iqtree-omp -nt 3 -st AA -s  rps12_aa.phy -m cpREV        -b 100 -pre rps12_aa_ts
iqtree-omp -nt 3 -st AA -s  rps14_aa.phy -m cpREV        -b 100 -pre rps14_aa_ts
iqtree-omp -nt 3 -st AA -s  rps15_aa.phy -m cpREV        -b 100 -pre rps15_aa_ts
iqtree-omp -nt 3 -st AA -s  rps16_aa.phy -m cpREV        -b 100 -pre rps16_aa_ts
iqtree-omp -nt 3 -st AA -s  rps18_aa.phy -m cpREV+G4     -b 100 -pre rps18_aa_ts
iqtree-omp -nt 3 -st AA -s  rps19_aa.phy -m cpREV        -b 100 -pre rps19_aa_ts
iqtree-omp -nt 3 -st AA -s  rps2_aa.phy  -m cpREV+G4     -b 100 -pre rps2_aa_ts
iqtree-omp -nt 3 -st AA -s  rps3_aa.phy  -m cpREV+G4     -b 100 -pre rps3_aa_ts
iqtree-omp -nt 3 -st AA -s  rps4_aa.phy  -m cpREV+G4     -b 100 -pre rps4_aa_ts
iqtree-omp -nt 3 -st AA -s  rps7_aa.phy  -m cpREV        -b 100 -pre rps7_aa_ts
iqtree-omp -nt 3 -st AA -s  rps8_aa.phy  -m cpREV        -b 100 -pre rps8_aa_ts
iqtree-omp -nt 3 -st AA -s  ycf1_aa.phy  -m cpREV+F+I+G4 -b 100 -pre ycf1_aa_ts
iqtree-omp -nt 3 -st AA -s  ycf2_aa.phy  -m cpREV+G4     -b 100 -pre ycf2_aa_ts
iqtree-omp -nt 3 -st AA -s  ycf3_aa.phy  -m cpREV        -b 100 -pre ycf3_aa_ts
iqtree-omp -nt 3 -st AA -s  ycf4_aa.phy  -m cpREV        -b 100 -pre ycf4_aa_ts
```

***

### ASTRAL-II

For input trees we used the 78 gene trees (x_genetrees.tre files) and their respective 100 bootstrap replicates (x-bs-file files) inferred using IQ-TREE.
```
#gene alignment
java -jar /home1/02700/deisejpg/bin/Astral/astral.4.11.1.jar -i ge_genetrees.tre -b ge-bs-file -o sp_aabs.tre 2> run2_gebs.log

#exon alignment
java -jar /home1/02700/deisejpg/bin/Astral/astral.4.11.1.jar -i ex_genetrees.tre -b ex-bs-file -o sp_aabs.tre 2> run2_exbs.log

#codon alignment
java -jar /home1/02700/deisejpg/bin/Astral/astral.4.11.1.jar -i cd_genetrees.tre -b cd-bs-file -o sp_aabs.tre 2> run2_cdbs.log

#amino acid alignment
java -jar /home1/02700/deisejpg/bin/Astral/astral.4.11.1.jar -i aa_genetrees.tre -b aa-bs-file -o sp_aabs.tre 2> run2_aabs.log
```
***

### SVDquartets

Nexus block for gene alignments:
```
#NEXUS

begin data;
    dimensions ntax=94 nchar=112750;
    format datatype=dna missing=? gap=- interleave=yes;
    matrix
    atpA_alignment ...
    atpB_alignment ...
    ...
;
end;

begin sets;
    charpartition loci =accD:1-4057,
                atpA:4058-5596,
                atpB:5597-7099,
                ...
;
end;

begin sets;
    taxpartition rosidsspecies =
        TAXON_1     : 1-1,
        TAXON_2     : 2-2,
        ...
;
end;

Begin paup;
Log File='/scratch/02700/deisejpg/1_chapter/svdquartets/Ge/concatenated_genes.log' start;
Outgroup 1 6 83;
SVDQuartets evalQuartets=all nthreads=48;
SaveTrees file='/scratch/02700/deisejpg/1_chapter/svdquartets/Ge/nobs.tre';
SVDQuartets bootstrap=standard treeFile='/scratch/02700/deisejpg/1_chapter/svdquartets/Ge/svdq_boot.tre'
;
end;
```

The same code was used for the exon alignment:
```
Begin paup;
Log File='/scratch/02700/deisejpg/1_chapter/svdquartets/Ex/concatenated_exon.log' start;
Outgroup 1 6 83;
SVDQuartets evalQuartets=all nthreads=48;
SaveTrees file='/scratch/02700/deisejpg/1_chapter/svdquartets/Ex/nobs_ex.tre';
SVDQuartets bootstrap=standard treeFile='/scratch/02700/deisejpg/1_chapter/svdquartets/Ex/svdq_boot_ex.tre'
;
end;
```

And for the codon alignment:
```
Begin paup;
Log File='/scratch/02700/deisejpg/1_chapter/svdquartets/Cd/concatenated_codon.log' start;
Outgroup 1 6 83;
SVDQuartets evalQuartets=all nthreads=48;
SaveTrees file='/scratch/02700/deisejpg/1_chapter/svdquartets/Cd/nobs_codon.tre';
SVDQuartets bootstrap=standard treeFile='/scratch/02700/deisejpg/1_chapter/svdquartets/Cd/svdq_cd_boot.tre'
;
end;
```
***
