dmaw

<b>Installation</b>: <b>Description</b>: dmaw is a software suite on the computation and application of a new distance based on absent words.

Pre-installation Instructions
=============================

   This installation requires the pre-installation of the cmake tool,
a C++ compiler, and the libraries: libdivsufsort and sdsl.

For Linux, you can install libraries libdivsufsort and sdsl via

 $ ./pre-install.sh 

Basic Instructions
==================

   The shell command `make -f Makefile.32-bit.gcc' should compile this 
program for 32-bit integers. 
   
   The shell command `make -f Makefile.64-bit.gcc' should compile this 
program for 64-bit integers. This requires double the amount of memory.

After compilation the binary `dmaw' will be created in the working 
directory, e.g. you may call it from this directory via

 $ ./dmaw -a DNA -i ./data/1.fasta -o 1.fasta.out

```
 Usage: dmaw
 Standard (Mandatory):
  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT'
                                      for protein  sequences. 
  -i, --input-file          <str>     (Multi)FASTA input filename.
  -o, --output-file         <str>     Output filename.
```

<b>Citations</b>:

```
When publishing work that is based on the results from dmaw please cite:

Castiglione G, Mantaci S, Restivo A. 
Some Investigations on Similarity Measures Based on Absent Words[J]. 
Fundamenta Informaticae, 2020, 171(1-4): 97-112.

```
<b>License</b>: GNU GPLv3 License; Copyright (C) 2021 Jia Gao
