
***
## I. GENOME ASSEMBLY
**Authors:** Deise JP Goncalves et al.  
**Date:** 2018/06/04  
**Citation:**   
### Goncalves, DJP, Simpson, BB, Ortiz, EM, Shimizu, GH, Jansen, RK. Incongruence between species tree and gene trees and phylogenetic signal variation in plastid genes. Molecular Phylogenetics and Evolution. In review
### DJP Goncalves, BB Simpson, GH Shimizu, RK Jansen, EM Ortiz. Genome assembly and phylogenomic data analyses using plastid data: contrasting species tree estimation methods. Data in Brief. Submitted
***

The pipeline for preprocessing the data, from raw data to assembled genomes, can be organized in directories with the following names:
```bash
1_raw
2_bbduk_clean_ar
3_bbduk_clean_phix
4_bbmerge_qtrim
5_ray_assembly
```

### `1_raw` 
In the current directory we store raw sequence data from genome skimming. Here we use _Oxalis drummondii_ as an example.
```bash
OxaliDG722_S34_L007_R1_001.fastq.gz
OxaliDG722_S34_L007_R2_001.fastq.gz
```

The commands bellow are prepared to run in each directory, which are represented by `2_bbduk_clean_ar`, `3_bbduk_clean_phix`, `4_bbmerge_qtrim`, and `5_ray_assembly`. The outputs are also written to be saved in the next directory, continuing the pipeline.

### `2_bbduk_clean_ar`
Below are the commands used for adaptor removal in bbduk from the BBTools suite.  
Installation guide at: [https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/].  
More information about `bbduk.sh` commands at: [https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/]
```bash
bbduk.sh in=../1_raw/OxaliDG722_S34_L007_R1_001.fastq.gz in2=../1_raw/OxaliDG722_S34_L007_R2_001.fastq.gz out=OxaliDG722_PE.ar.fq.gz outs=OxaliDG722_SE.ar.fq.gz ktrim=r k=23 hdist=1 mink=11 hdist2=1 tbo stats=OxaliDG722.stats.txt ref=$HOME/bin/bbmap/resources/adapters.fa &> OxaliDG722_bbduk_ar.log
```

### `3_bbduk_clean_phix`
This step is also done in `bbduk.sh` and it is done to remove PhiX matches from the sequences.
```bash
bbduk.sh in=../2_bbduk_clean_ar/OxaliDG722_PE.ar.fq.gz out=OxaliDG722_PE.cl.fq.gz outs=OxaliDG722_SE.cl.fq.gz k=31 hdist=1 stats=OxaliDG722.statsPhix.txt ref=$HOME/bin/bbmap/resources/phix174_ill.ref.fa.gz &> OxaliDG722.bbduk_cl.log
```

### `4_bbmerge_qtrim`
`bbmerge.sh` is another software from the `BBTools` suite and is used to overlap paired reads into single reads. 
More information about the commands at: [https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/]
```bash
bbmerge.sh in=../3_bbduk_clean_phix/OxaliDG722_PE.cl.fq.gz out=OxaliDG722_ME.cl.fq.gz outu=OxaliDG722_PE.cl.me.fq.gz qtrim2=t trimq=10 maq=10 ml=61 &> OxaliDG722.bbmerge.log
```

### `5_ray_assembly`
In the current directory you must have the references that will be used by the assembler to separate contigs that maps to each genomic compartment. As called in Ray code below, "-search refs/ptd", the "refs" directory should be saved within the current directory and should contain three other directories, "mit", "ptd", and "rib"; standing for mitochondria, plastid, and ribosome, respectively.
```
5_ray_assembly+
              |-refs+
                    |-mit+
                         |-mit.fasta
                    |-ptd+
                         |-ptd.fasta
                    |-rib+
                         |-rib.fasta
```
The reference files were prepared using reference data from species available on NCBI.  
Below is an example of the command to run `Ray` with k-mer 55 to k-mer 135.  
More information about how to download and about the commands at [http://denovoassembler.sourceforge.net/index.html]
```bash
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 55 -o OxaliDG722_PE_K55_RAY &>OxaliDG722_PE_K55_RAY.log
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 65 -o OxaliDG722_PE_K65_RAY &>OxaliDG722_PE_K65_RAY.log
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 75 -o OxaliDG722_PE_K75_RAY &>OxaliDG722_PE_K75_RAY.log
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 85 -o OxaliDG722_PE_K85_RAY &>OxaliDG722_PE_K85_RAY.log
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 95 -o OxaliDG722_PE_K95_RAY &>OxaliDG722_PE_K95_RAY.log
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 105 -o OxaliDG722_PE_K105_RAY &>OxaliDG722_PE_K105_RAY.log
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 115 -o OxaliDG722_PE_K115_RAY &>OxaliDG722_PE_K115_RAY.log
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 125 -o OxaliDG722_PE_K125_RAY &>OxaliDG722_PE_K125_RAY.log
ibrun Ray -i ../4_bbmerge_qtrim/OxaliDG722_PE.cl.me.fq.gz -s ../4_bbmerge_qtrim/OxaliDG722_ME.cl.fq.gz -search refs/ptd -search refs/mit -search refs/rib -k 135 -o OxaliDG722_PE_K135_RAY &>OxaliDG722_PE_K135_RAY.log
```

After preprocessing the data, output files were downloaded from the cluster. Plastomes were completed and sequences of all plastid protein-coding genes were extracted using Geneious v. 7.1.9 (Kearse et al., 2012).
