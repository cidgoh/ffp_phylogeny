Feature Frequency Profile Phylogenies
=====================================


Introduction
------------

FFP (Feature frequency profile) is an alignment free comparison tool for phylogenetic analysis and text comparison. It can be applied to nucleotide sequences, complete genomes, proteomes and even used for text comparison.  This software is a Galaxy (http://galaxyproject.org) tool for calculating FFP on one or more fasta sequence or text datasets.

The original command line ffp-phylogeny code is at http://ffp-phylogeny.sourceforge.net/ .  This tool uses Aaron Petkau's modified version: https://github.com/apetkau/ffp-3.19-custom .  Aaron has quite a good writeup of the technique as well at https://github.com/apetkau/microbial-informatics-2014/tree/master/labs/ffp-phylogeny .

**Installation Note** : Your Galaxy server will need the groff package to be installed on it first (to generate ffp-phylogeny man pages).  A cryptic error will occur if it isn't: "troff: fatal error: can't find macro file s".  This is different from the "groff-base" package.

This Galaxy tool prepares a mini-pipeline consisting of **[ffpry | ffpaa | ffptxt] > [ ffpfilt | ffpcol > ffprwn] > ffpjsd > ffptree**  .  The last step is optional - by deselecting the "Generate Tree Phylogeny" checkbox, the tool will output a distance matrix rather than a Newick (.nhx) formatted tree file.

Each sequence or text file has a profile containing tallies of each feature found.  A feature is a string of valid characters of given length. 

For nucleotide data, by default each character (ATGC) is grouped as either purine(R) or pyrmidine(Y) before being counted.

For amino acid data, by default each character is grouped into one of the following: (ST),(DE),(KQR),(IVLM),(FWY),C,G,A,N,H,P. Each group is represented by the first character in its series.

One other key concept is that a given feature, e.g. "TAA" is counted in forward AND reverse directions, mirroring the idea that a feature's orientation is not so important to distinguish when it comes to alignment-free comparison.  The counts for "TAA" and "AAT" are merged.
 
The labeling of the resulting counted feature items is perhaps the trickiest concept to master.  Due to computational efficiency measures taken by the developers, a feature that we see on paper as "TAC" may be stored and labeled internally as "GTA", its reverse compliment.  One must look for the alternative if one does not find the original. 

Also note that in amino acid sequences the stop codon "*" (or any other character that is not in the Amino acid alphabet) causes that character frame not to be counted.  Also, character frames never span across fasta entries.

A few tutorials:
 * http://sourceforge.net/projects/ffp-phylogeny/files/Documentation/tutorial.pdf
 * https://github.com/apetkau/microbial-informatics-2014/tree/master/labs/ffp-phylogeny

-------
**Note**

Taxonomy label details: If each file contains one profile, the file's name is used to label the profile.  If each file contains fasta sequences to profile individually, their fasta identifiers will be used to label them.  The "short labels" option will find the shortest label that uniquely identifies each profile.  Either way, there are some quirks: ffpjsd clips labels to 10 characters if they are greater than 50 characters, so all labels are trimmed to 50 characters first.  Also "id" is prefixed to any numeric label since some tree visualizers won't show purely numeric labels.  In the accidental case where a Fasta sequence label is a duplicate of a previous one it will be prefixed by "DupLabel-".

The command line ffpjsd can hang if one provides an l-mer length greater than the length of file content.  One must identify its process id ("ps aux | grep ffpjsd") and kill it ("kill [process id]").

Finally, it is possible for the ffptree program to generate a tree where some of the branch distances are negative. See https://www.biostars.org/p/45597/

-------
**References**
 
The development of the ffp-phylogeny command line software should be attributed to:

Sims GE, Jun S-R, Wu GA, Kim S-H. Alignment-free genome comparison with feature frequency profiles (FFP) and optimal resolutions. Proceedings of the National Academy of Sciences of the United States of America 2009;106(8):2677-2682. doi:10.1073/pnas.0813249106.

