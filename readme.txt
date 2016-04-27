Notice:
This is the source code fold.
The largest number of contigs are 2^31. The input file must be FASTA format. Please read Requirement and How to use this tool before running.

 Thanks for reporting bugs and unexpected output.

Requriement:

This source code is suitable for all unix-like 64-bit system with gcc installed.

How to use this method:
(1)make  (the executive file 'IFCM' will be in this dictionary)
(2)use Nonpareil (https://github.com/lmrodriguezr/nonpareil/) or other similar method to estimate the coverage of metagenome
(3)./IFCM inputfile coverage
	for example: ./IFCM example.fna coverage

Output:
result.txt file with the cluster label of each contig, from Cluster 0 to Cluster k.

Contact us:
laoniu313@qq.com