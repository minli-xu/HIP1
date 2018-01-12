Two data files are hosted.

[trans.1368.txt] has the 1368 transcripts from Vijayan et al. 2011.
From the original 1415 transcripts. 47 shorter transcripts that are
completely contained in longer transcripts on the same strand, are 
removed, resulting 1368 transcripts.

[trans.1149.txt] has the 1149 transcripts from Vijayan et al. 2011.
From the set of 1368 transcripts. We attached additional information
about the transcripts, such as ACEu and CAI. If these additional
information could not be linked to a transcript in the 1368 set,
we excluded these transcripts, resulting 1149 transcripts. 


Column information for [trans.1368.txt]

Column 1	Transcript Start
Column 2	Transcript End
Column 3	Strand. 1 for positive(+) strand and 0 for negative(-) strand
Column 4	Number of ORF contained
Column 5	5'- untranslated region (UTR) length
Column 6	3'- untranslated region (UTR) length


Column information for [trans1149.txt]

Column 1	Transcript Start
Column 2	Transcript End
Column 3	Strand. 1 for positive(+) strand and 0 for negative(-) strand
Column 4	Number of ORF contained
Column 5	5'- untranslated region (UTR) length
Column 6	3'- untranslated region (UTR) length
Column 7*	Mean RNA sequencing over transcript (w1)
Column 8**	Molecules per cell (mRNA/cell) (w2)
Column 9	Average protein length (unit:aa)
Column 10	ORF average GC content
Column 11	ACEu
Column 12	ACEz
Column 13	CAI
Column 14	CBI	

* w1.  "mean of the raw RNA sequencing reads over the full transcript"
**w2.  "number of transcripts per cell assuming a total of 1,500 mRNAs per cell"


Note: the transcript's coordinates, as well as transcript abundance (w1 and w2)
are from this study:

Vikram Vijayan, Isha H Jain, Erin K O’Shea
A high resolution map of a cyanobacterial transcriptome
Genome Biology, 2011, Volume 12, Number 5, Page 1

https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r47