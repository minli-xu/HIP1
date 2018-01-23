

The periodicity of the spatial motif distribution was assessed by
fitting a damped sine wave to the autocorrelation of the histogram
of inter-motif distances (spacings), following a procedure described
in ... (URL OF THE SUPPLEMENTARY TEXT TO BE FILLED LATER, WITH A 
SPECIFIC PAGE OR SECTION NUMBER)

The results of the fitting procedure depends on the bin size used to
aggregate the inter-motif spacings. The fitting procedure was carried
out for 280 bin sizes ranging from 25bp to 7000bp in 25bp intervals.

For each genome, this procedure was carried out (1) for the genuine
data, based on the coordinates HIP1-like (HIP1 + HIP1*) motifs and
(2) for 200 "randomized" genomes, wherein motif sites are randomly
reassigned within the genome, with the restriction that sites must
be separated by an interval that is not shorter than the length of
the motif.

The parameters of the damped sine wave are 

 - Goodness of fit (GoF)
 - Period, 
 - Amplitude


======================  HIP1-like motifs ======================

The directory HIP1_motifs_compilation_aug3/ contains the parameters of
the inferred periodic signal for HIP1-like motifs in 50 HIP1 rich genomes. 

There are six files for each genome, three for the genuine motif coordinates
and three for the randomized coordinates.

The file names for the genuine coordinates are of the form

  genome.motif_type.parameter_type.

These files contain a single row of parameter values, where the i-th
entry in the row corresponds to the best-fit sine wave obtained with
bin of size i*25 bp.

For example, the file Anabaena_90_uid179383.HIP1.GoF contains one
row of 280 real numbers (space delimited), corresponding to the
goodness-of-fit values obtained with bins of size 25bp, 50bp,
... 7000bp.  

The file names for the randomized coordinates are of the form

  genome.motif_type.parameter_type.random.

Each of these files contains 200 rows, one row for each randomized
genome.  The format for the individual rows is the same as for the
genuine data.



======================  Control motifs ======================

The directory Control_motifs_compilation_aug3/ contains the parameters of
the inferred periodic signal for control motifs in 46 HIP1 rich genomes. 

The files are named and formatted similarly as the files for HIP1-like
motifs in HIP1_motifs_compilation_aug3/.

* Note that there are totally six Synechocystis sp. PCC6803 genomes in our study.
Only the control motifs from the genome of Synechocystis_PCC_6803_uid57659 
was used for the analysis, as the six Synechocystis sp. PCC6803 genomes are
extremely similar.
