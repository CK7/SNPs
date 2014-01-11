SNPs
====

This is the code for the program SNPs that computes base calls information from one or more SAM files for sequences in the provided fasta file.

Installation
------------

Simply run "make" under the directory to which the SNPs files were downloaded.
The program should now be ready to run under bin/SNPs.

Running SNPs
------------
You can get the following information by running SNPs without parameters:<br>
```
$ bin/SNPs

SNPs version 1.02, 12/31/13

Computes base calls, insertions and deletions for all positions on every sequence in the input files.

Usage: bin/SNPs -assembly <sequence-file> -out <output-file> [-snp <snp-threshold>] [--silent]
                [-sam  <read-mapping-file1> <read-mapping-file2> ... <read-mapping-filen>]
Where
  -assembly         :  path to assembly file
  -out              :  path for the output file
  -sam              :  mapping files are in SAM format (only option in current version)
  -snp              :  threshold for determining SNPs (default: 0.2)
  --silent          :  progrsm will do its work without status messages to stderr (default: verbose)

OUTPUT
  output consists of two type of lines:
  - Header lines    :  mark the beginning of a new sequence. Format is
                    :  ><reference-seq-name> /size=<size> /coverage=<coverage> /gc=<%G+C> /Ns=<# of Ns>
  - Base call lines :  present information for a certain position of the sequence whose header was last seen.
                    :  Format is
                    :  <pos> <reference-char> <coverage> <# A's> <# C's> <# G's> <# T's> <# N's> <# Insertions> <# Deletions> <additional-columns>
                    :  where <additional-columns> are either empty or contains the following:
                    :  + X -> Y, if X is the reference char but Y is the consensus
                    :  + SNP if the number of calls of the non-consensus char is more than <SNP-threshold> of the A, C, G and T (but not N) calls
                    :  + INSERTION if the number of insertion calls is more than <SNP-threshold> of A, C, G, T AND N calls for the position
                    :  + DELETION if the number of deletion calls if more than <SNP-threshold> of A, C, G, T, N AND deletion call for the position
                    :  All fields, including additional columns, are separated by tabs.
                    :  Insertions refer to insertions coming right after the position in which they are recognized. Insertions can be of 1 or more
                    :  bases, in this version only the presence of an insertion is reported.

Please report bugs to itai.sharon@gmail.com

```
SNPs should be able to take any legal SAM file. It was tested on output from bowtie and bowtie2.

Output file looks as follows:
```
Header lines:
><reference-seq-name> /size=<size> /coverage=<coverage> /gc=<%G+C> /Ns=<# of Ns>

Base call lines:
<pos> <reference-char> <coverage> <# A's> <# C's> <# G's> <# T's> <# N's> <# Insertions> <# Deletions> <additional-columns>
```
See above description for the different columns.

For example:
```
>contig1 /size=197362 /coverage=853.933 /gc=0.343197 /Ns=0
1	C	435	29	394	4	8	0	3	0				
2	G	443	15	8	406	14	0	4	0				
3	C	447	4	409	14	20	0	10	0				
4	A	450	412	9	13	16	0	38	0				
5	A	459	454	1	3	1	0	4	0				
6	A	461	455	4	1	1	0	0	0				
7	C	461	1	458	1	1	0	2	0				
8	C	467	1	465	0	1	0	0	0				
9	C	472	0	471	0	1	0	0	0
:
2700	A	859	6	846	2	5	0	2	0	A -> C
:
146733	G	750	30	153	559	8	0	0	0		SNP
:
146745	C	737	498	231	4	4	0	1	1	C -> A	SNP
:
160293	A	894	879	6	0	9	0	840	0			INSERTION
:
```
