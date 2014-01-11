SNPs
====

This is the code for the program SNPs that generates a report for the coverage and the occurences of each DNA character in every position on the input sequence(s).

Installation
------------
1. Download all files in the general_purpose_dna_cpp_classes reporsitory and store them under a directory named "common" (otherwise you'll have to modify the include 
directories in this repository).
2. Download the files in this directory and store them under directory "SNPs".
3. Create a directory called "bin" under "SNPs".
4. run the command 
   make -f makefile.SNPs

The program should now be ready to run under bin/SNPs.

Execution
---------
Command line looks as follows:<br>
```
$ bin/SNPs
Usage: bin/SNPs -assembly <assembly file> -out <output-file> [--phred33-quals|--phred64-quals] [--silent]
                [-sam  \<read-mapping-file1> <read-mapping-file2> ... <read-mapping-filen>]
Where
  -assembly       :  path to assembly file
  -out            :  path for the output file
  -sam            :  mapping files are in SAM format (only option in current version)
  --phred33-quals :  Quality scores of reads are ASCII characters equal to the Phred quality plus 33 (same as bowtie)
  --phred64-quals :  Quality scores of reads are ASCII characters equal to the Phred quality plus 64 (default, same as bowtie)
  --silent        :  progrsm will do its work without status messages to stderr (default: verbose)
```
I use bowtie for generating the mapping files.

Output file looks as follows:
```
<Sequence-name>   <position>  <scaf-character>  <coverage>  <#A's>  <#C's>  <#G's>  <#T's>  <#N's>  <replacement> <SNP>
```
Where
- replacement is reported if the most abundant char is not the one reported on the scaffold
- SNP is reported if the most abundant char is < 80% of the chars in that position
- Both are considered only for positions with coverage >= 10.

e.g.
```
>NODE_1_length_197261_cov_539.935364 /size=197362 /coverage=802.168 /gc=0.343278 /Ns=0
0  C	375	4	370	0	0	1
1	G	378	1	0	377	0	0
2	C	382	0	382	0	0	0
3	A	383	383	0	0	0	0
4	A	392	392	0	0	0	0
5	A	394	392	0	1	1	0
  :
4946  C	65	11	50	1	3	0		SNP
4947  T	67	14	7	7	39	0		SNP
  :
9504  G	920	609	60	92	158	1	G -> A	SNP
  :
```

fix_SNPs.pl
-----------

This script takes a fasta file as well as a SNPs report for the sequences in the file and changes all positions in which the most abundant char is different from the one 
appearing in the sequence.
Usage:
```
$ ./fix_SNPs.pl 

Usage: ./fix_snps.pl <fasta-file> <SNPs report> <fasta-file-out>

```
Output is written to fasta-file-out. All changed positions will be written to stderr.
