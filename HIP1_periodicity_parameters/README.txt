Parameters estimated from sine wave fitting
------------------------------------------------------

The directiory HIP1_motifs_compilation_aug3/ has the Goodness of fit, inferred 
Period, and Amplitude parameters, for HIP1-like(HIP1 + HIP1*) data

The file name is in the form of:  genome.motif_type.parameter_type,
for example:

Anabaena_90_uid179383.HIP1.GoF

or like this for random data

Anabaena_90_uid179383.HIP1.GoF.random

Each row is from one of the 200 randomization of the coordinates. So each random 
file has 200 rows, and 280 columns for 280 bins (25 to 7000 bp with 25 bp 
increasement). 

Each line is space delimited.

The non-random one is similar. 1 row and 280 bins.

-----------------------------------------------------------------------------

The directory Control_motifs_compilation_aug3/ is similar, but for control motifs.

The file name is in the form of : genome.motif_type.motif_type.parameter_type,
for example:

Anabaena_90_uid179383.AGGCCT.Amp

or like this for random data

Anabaena_90_uid179383.AGGCCT.Amp.random 

The format of the file is same as for HIP1 motifs (200x280 matrix for random, and
1 x 280 matrix for non-random)

Each line is space delimited.


* Note that there are totally six Synechocystis sp. PCC6803 genomes in our study.
Only the hexamers from the genome of Synechocystis_PCC_6803_uid57659 was used 
for the analysis, as the six Synechocystis sp. PCC6803 genomes are extremely 
similar.
