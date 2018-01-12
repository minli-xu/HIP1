# 9/8/2015 11:53:25 AM added the alignment score (id% and coverage)
# 9/15/2015 2:08:45 PM be able to do mismatch 1  . perl compute_consv_SAK.pl Two_Thermo_BP1_NK55a_gai GCGATCGC 1

use strict;
use warnings;

my $input_aln = shift; # the MAUVE algnment file
#my $input_motif = 'GCGATCGC'; # can be ys motifs or control.
my $input_motif = shift; # can be ys motifs or control.
my $flagm1 = shift;
$flagm1 = 0 unless defined $flagm1;


if ((!defined $input_motif)||($input_motif eq '')){
	print "PLEASE INPUT THE INPUT MOTIF !!! \n";
}elsif (($input_motif eq 'c') || ($input_motif eq 'C')){
	
	my @all_c = ('AGGGCCCT','GAGGCCTC','GGAGCTCC','GGGATCCC','AGGCGCCT','GAGCGCTC','GGACGTCC','GGCATGCC','AGCGCGCT','GACGCGTC','GCAGCTGC','AGCCGGCT','GACCGGTC','GCACGTGC','GCCATGGC','ACGGCCGT','CAGGCCTG','CGAGCTCG','CGGATCCG','ACGCGCGT','CAGCGCTG','CGACGTCG','CGCATGCG','ACCGCGGT','CACGCGTG','CCAGCTGG','CCGATCGG','ACCCGGGT','CACCGGTG','CCACGTGG','CCCATGGG','TGGGCCCA','GTGGCCAC','GGTGCACC','GGGTACCC','TGGCGCCA','GTGCGCAC','GGTCGACC','GGCTAGCC','TGCGCGCA','GTCGCGAC','GCTGCAGC','GCGTACGC','TGCCGGCA','GTCCGGAC','GCTCGAGC','GCCTAGGC','TCGGCCGA','CTGGCCAG','CGTGCACG','CGGTACCG','TCGCGCGA','CTGCGCAG','CGTCGACG','CGCTAGCG','TCCGCGGA','CTCGCGAG','CCTGCAGG','CCGTACGG','TCCCGGGA','CTCCGGAG','CCTCGAGG','CCCTAGGG');
	if ($input_aln eq 'syh_syg_gai'){
		@all_c = ('AGGGCCCT','GAGGCCTC','GGAGCTCC','GCGATCGC','AGGCGCCT','GAGCGCTC','GGACGTCC','GGCATGCC','AGCGCGCT','GACGCGTC','GCAGCTGC','AGCCGGCT','GACCGGTC','GCACGTGC','GCCATGGC','ACGGCCGT','CAGGCCTG','CGAGCTCG','CGGATCCG','ACGCGCGT','CAGCGCTG','CGACGTCG','CGCATGCG','ACCGCGGT','CACGCGTG','CCAGCTGG','CCGATCGG','ACCCGGGT','CACCGGTG','CCACGTGG','CCCATGGG','TGGGCCCA','GTGGCCAC','GGTGCACC','GGGTACCC','TGGCGCCA','GTGCGCAC','GGTCGACC','GGCTAGCC','TGCGCGCA','GTCGCGAC','GCTGCAGC','GCGTACGC','TGCCGGCA','GTCCGGAC','GCTCGAGC','GCCTAGGC','TCGGCCGA','CTGGCCAG','CGTGCACG','CGGTACCG','TCGCGCGA','CTGCGCAG','CGTCGACG','CGCTAGCG','TCCGCGGA','CTCGCGAG','CCTGCAGG','CCGTACGG','TCCCGGGA','CTCCGGAG','CCTCGAGG','CCCTAGGG');
	}
	my $all_found_in_ref = 0;
	my $all_found_in_tar = 0;
	my $all_consv_in_tar = 0;
	foreach my $ctrlmotif (@all_c){
		my ($found_in_ref,$found_in_tar,$consv_in_tar) = compute_consv($ctrlmotif);
		$all_found_in_ref += $found_in_ref;
		$all_found_in_tar += $found_in_tar;
		$all_consv_in_tar += $consv_in_tar;	
		print "$ctrlmotif\t$found_in_ref\t$found_in_tar\t$consv_in_tar\n";
	}
	#print "\nCONTROL\t$all_found_in_ref\t$all_found_in_ref\t$all_consv_in_tar\n";	

	
	my ($found_in_ref,$found_in_tar,$consv_in_tar) = compute_consv('GCGATCGC');	
	my $C_score = sprintf("%.3f",$consv_in_tar / $found_in_ref);
	my $S_score = sprintf("%.3f",$consv_in_tar / ($found_in_ref+$found_in_tar-$consv_in_tar));
	print "\nGCGATCGC\t$found_in_ref\t$found_in_tar\t$consv_in_tar\t$C_score\t$S_score\n";
	
	my $c_C_score = sprintf("%.3f",$all_consv_in_tar / $all_found_in_ref);
	my $c_S_score = sprintf("%.3f",$all_consv_in_tar / ($all_found_in_ref+$all_found_in_tar-$all_consv_in_tar));
	print "CONTROL\t$all_found_in_ref\t$all_found_in_tar\t$all_consv_in_tar\t$c_C_score\t$c_S_score\n";

}elsif (($input_motif eq 'm1c') || ($input_motif eq 'm1C')){
	
	my @all_c = ('AGGGCCCT','GAGGCCTC','GGAGCTCC','GGGATCCC','AGGCGCCT','GAGCGCTC','GGACGTCC','GGCATGCC','AGCGCGCT','GACGCGTC','GCAGCTGC','AGCCGGCT','GACCGGTC','GCACGTGC','GCCATGGC','ACGGCCGT','CAGGCCTG','CGAGCTCG','CGGATCCG','ACGCGCGT','CAGCGCTG','CGACGTCG','CGCATGCG','ACCGCGGT','CACGCGTG','CCAGCTGG','CCGATCGG','ACCCGGGT','CACCGGTG','CCACGTGG','CCCATGGG','TGGGCCCA','GTGGCCAC','GGTGCACC','GGGTACCC','TGGCGCCA','GTGCGCAC','GGTCGACC','GGCTAGCC','TGCGCGCA','GTCGCGAC','GCTGCAGC','GCGTACGC','TGCCGGCA','GTCCGGAC','GCTCGAGC','GCCTAGGC','TCGGCCGA','CTGGCCAG','CGTGCACG','CGGTACCG','TCGCGCGA','CTGCGCAG','CGTCGACG','CGCTAGCG','TCCGCGGA','CTCGCGAG','CCTGCAGG','CCGTACGG','TCCCGGGA','CTCCGGAG','CCTCGAGG','CCCTAGGG');
	if ($input_aln eq 'syh_syg_gai'){
		@all_c = ('AGGGCCCT','GAGGCCTC','GGAGCTCC','GCGATCGC','AGGCGCCT','GAGCGCTC','GGACGTCC','GGCATGCC','AGCGCGCT','GACGCGTC','GCAGCTGC','AGCCGGCT','GACCGGTC','GCACGTGC','GCCATGGC','ACGGCCGT','CAGGCCTG','CGAGCTCG','CGGATCCG','ACGCGCGT','CAGCGCTG','CGACGTCG','CGCATGCG','ACCGCGGT','CACGCGTG','CCAGCTGG','CCGATCGG','ACCCGGGT','CACCGGTG','CCACGTGG','CCCATGGG','TGGGCCCA','GTGGCCAC','GGTGCACC','GGGTACCC','TGGCGCCA','GTGCGCAC','GGTCGACC','GGCTAGCC','TGCGCGCA','GTCGCGAC','GCTGCAGC','GCGTACGC','TGCCGGCA','GTCCGGAC','GCTCGAGC','GCCTAGGC','TCGGCCGA','CTGGCCAG','CGTGCACG','CGGTACCG','TCGCGCGA','CTGCGCAG','CGTCGACG','CGCTAGCG','TCCGCGGA','CTCGCGAG','CCTGCAGG','CCGTACGG','TCCCGGGA','CTCCGGAG','CCTCGAGG','CCCTAGGG');
	}
	
	my $all_found_in_ref = 0;
	my $all_found_in_tar = 0;
	my $all_consv_two_m1_S = 0;
	my $all_consv_two_m1_NS = 0;
	my $all_consv_m1_m0 = 0;
	my $all_notconsv_m1_x = 0;
	my $all_consv_two_m0 = 0;
	
	foreach my $ctrlmotif (@all_c){
		my ($found_in_ref,$found_in_tar,$consv_two_m1_S,$consv_two_m1_NS,$consv_m1_m0,$notconsv_m1_x,$consv_two_m0);
	   	   ($found_in_ref,$found_in_tar,$consv_two_m1_S,$consv_two_m1_NS,$consv_m1_m0,$notconsv_m1_x,$consv_two_m0) = compute_consv_3($ctrlmotif);
		#my $C_score = sprintf("%.3f",$consv_in_tar / $found_in_ref);
		#my $S_score = sprintf("%.3f",$consv_in_tar / ($found_in_ref+$found_in_tar-$consv_in_tar));
		print "$input_aln\t$input_motif\t$found_in_ref\t$found_in_tar\t$consv_two_m1_S\t$consv_two_m1_NS\t$consv_m1_m0\t$notconsv_m1_x\t$consv_two_m0\n";	
		$all_found_in_ref    +=  $found_in_ref   ;
		$all_found_in_tar    +=  $found_in_tar   ;
		$all_consv_two_m1_S  +=  $consv_two_m1_S ;
		$all_consv_two_m1_NS +=  $consv_two_m1_NS;
		$all_consv_m1_m0     +=  $consv_m1_m0    ;
		$all_notconsv_m1_x   +=  $notconsv_m1_x  ;
		$all_consv_two_m0    +=  $consv_two_m0   ;
	
	}
	
	print "\nCONTROL\t$all_found_in_ref\t$all_found_in_tar\t$all_consv_two_m1_S\t$all_consv_two_m1_NS\t$all_consv_m1_m0\t$all_notconsv_m1_x\t$all_consv_two_m0\n";

}else{
	if ($flagm1 == 0){
		my ($found_in_ref,$found_in_tar,$consv_in_tar,$n_pairs,$idperc_match,$idperc_total);
		   ($found_in_ref,$found_in_tar,$consv_in_tar,$n_pairs,$idperc_match,$idperc_total) = compute_consv_2($input_motif);;
		my $C_score = sprintf("%.3f",$consv_in_tar / $found_in_ref);
		my $S_score = sprintf("%.3f",$consv_in_tar / ($found_in_ref+$found_in_tar-$consv_in_tar));
		print "$input_motif\t$found_in_ref\t$found_in_tar\t$consv_in_tar\t$C_score\t$S_score\t$n_pairs\t$idperc_match\t$idperc_total\n";
	}elsif ($flagm1 == 1){
		my ($found_in_ref,$found_in_tar,$consv_two_m1_S,$consv_two_m1_NS,$consv_m1_m0,$notconsv_m1_x,$consv_two_m0);
		   ($found_in_ref,$found_in_tar,$consv_two_m1_S,$consv_two_m1_NS,$consv_m1_m0,$notconsv_m1_x,$consv_two_m0) = compute_consv_3($input_motif);;
		#my $C_score = sprintf("%.3f",$consv_in_tar / $found_in_ref);
		#my $S_score = sprintf("%.3f",$consv_in_tar / ($found_in_ref+$found_in_tar-$consv_in_tar));
		print "$input_aln\t$input_motif\t$found_in_ref\t$found_in_tar\t$consv_two_m1_S\t$consv_two_m1_NS\t$consv_m1_m0\t$notconsv_m1_x\t$consv_two_m0\n";
	}
}


sub compute_consv_3{
	
	my $this_motif = shift;
	
	my $l = length $this_motif;
	
	my $head_line;
	my $big_seq_str;
	my @all_seqs;
	open(IN,$input_aln) or die("can not open aln file \n");
	while(<IN>){
		my $line =$_;
		chomp $line;
		next if ($line =~/^#/);
		next if ($line =~/^=/);
		if ($line =~/^>/){
			push @all_seqs,$big_seq_str;
			$big_seq_str = '';
		}else{
			$big_seq_str = $big_seq_str.$line;
		}
	}
	push @all_seqs,$big_seq_str;
	shift @all_seqs;
	close IN;
	#print "$#all_seqs +1 seqs obtained ...\n" unless ($input_motif eq 'c');
	
	my $n_pairs = scalar (@all_seqs) /2;
	
	my $found_in_ref = 0;
	my $found_in_tar = 0;
	my $consv_two_m1_S = 0;
	my $consv_two_m1_NS = 0;
	my $consv_m1_m0 = 0;
	my $consv_two_m0 = 0;
	my $notconsv_m1_x = 0;
	my $idperc_match = 0;
	my $idperc_total = 0;
	
	
	for my $i (0 .. ($n_pairs-1)){
		my $ref_g = $all_seqs[$i*2];
		my $tar_g = $all_seqs[$i*2+1];
		
		#my $found_in_ref = 0;
		#my $consv_in_tag = 0;
		
		my $n1 = length $ref_g;
		my $n2 = length $tar_g;
		
		die("not equal length in ref and tar seq !!!\n$n1 $n2 \n") unless (length $ref_g eq length $tar_g);
		my $L = length $ref_g;
		
		for my $j (0 .. ($L - $l)){
			my $substr_ref = substr($ref_g,$j,$l);
			my $substr_tar = substr($tar_g,$j,$l);
			if((mismatch($substr_ref,$this_motif)==1) && (mismatch($substr_tar,$this_motif)==1)){
				$found_in_ref++;
				$found_in_tar++;
				$consv_two_m1_S  ++ if ($substr_ref eq $substr_tar);
				$consv_two_m1_NS ++ if ($substr_ref ne $substr_tar);
			}elsif ((mismatch($substr_ref,$this_motif)==1) && (mismatch($substr_tar,$this_motif)==0)){
				$found_in_ref++;
				$consv_m1_m0 ++;	
			}elsif ((mismatch($substr_ref,$this_motif)==0) && (mismatch($substr_tar,$this_motif)==1)){
				$found_in_tar++;
				$consv_m1_m0 ++;	
			}elsif ((mismatch($substr_ref,$this_motif)==0) && (mismatch($substr_tar,$this_motif)==0)){
				$consv_two_m0 ++;
			}elsif ((mismatch($substr_ref,$this_motif)==1) && (mismatch($substr_tar,$this_motif)>=2)){
				$found_in_ref++;
				$notconsv_m1_x ++;	
			}elsif ((mismatch($substr_ref,$this_motif)>=2) && (mismatch($substr_tar,$this_motif)==1)){
				$found_in_tar++;
				$notconsv_m1_x ++;	
			}
			#$idperc_match ++ if (substr($ref_g,$j,1)) eq (substr($tar_g,$j,1));
			#$idperc_total ++;
		}
		
	}

	return ($found_in_ref,$found_in_tar,$consv_two_m1_S,$consv_two_m1_NS,$consv_m1_m0,$notconsv_m1_x,$consv_two_m0); 
}


sub compute_consv_2{
	
	my $this_motif = shift;
	
	my $l = length $this_motif;
	
	my $head_line;
	my $big_seq_str;
	my @all_seqs;
	open(IN,$input_aln) or die("can not open aln file \n");
	while(<IN>){
		my $line =$_;
		chomp $line;
		next if ($line =~/^#/);
		next if ($line =~/^=/);
		if ($line =~/^>/){
			push @all_seqs,$big_seq_str;
			$big_seq_str = '';
		}else{
			$big_seq_str = $big_seq_str.$line;
		}
	}
	push @all_seqs,$big_seq_str;
	shift @all_seqs;
	close IN;
	#print "$#all_seqs +1 seqs obtained ...\n" unless ($input_motif eq 'c');
	
	my $n_pairs = scalar (@all_seqs) /2;
	
	my $found_in_ref = 0;
	my $found_in_tar = 0;
	my $consv_in_tar = 0;
	my $idperc_match = 0;
	my $idperc_total = 0;
	
	
	for my $i (0 .. ($n_pairs-1)){
		my $ref_g = $all_seqs[$i*2];
		my $tar_g = $all_seqs[$i*2+1];
		
		#my $found_in_ref = 0;
		#my $consv_in_tag = 0;
		
		my $n1 = length $ref_g;
		my $n2 = length $tar_g;
		
		die("not equal length in ref and tar seq !!!\n$n1 $n2 \n") unless (length $ref_g eq length $tar_g);
		my $L = length $ref_g;
		
		for my $j (0 .. ($L - $l)){
			my $substr_ref = substr($ref_g,$j,$l);
			my $substr_tar = substr($tar_g,$j,$l);
			if(($substr_ref eq $this_motif) && ($substr_tar eq $this_motif)){
				$found_in_ref ++;
				$found_in_tar ++;
				$consv_in_tar ++;
			}elsif ($substr_ref eq $this_motif){
				$found_in_ref ++;
			}elsif ($substr_tar eq $this_motif){
				$found_in_tar ++;	
			}
			$idperc_match ++ if (substr($ref_g,$j,1)) eq (substr($tar_g,$j,1));
			$idperc_total ++;
		}
		
	}


	return ($found_in_ref,$found_in_tar,$consv_in_tar,$n_pairs,$idperc_match,$idperc_total); 
}


sub compute_consv{
	
	my $this_motif = shift;
	
	my $l = length $this_motif;
	
	my $head_line;
	my $big_seq_str;
	my @all_seqs;
	open(IN,$input_aln) or die("can not open aln file \n");
	while(<IN>){
		my $line =$_;
		chomp $line;
		next if ($line =~/^#/);
		next if ($line =~/^=/);
		if ($line =~/^>/){
			push @all_seqs,$big_seq_str;
			$big_seq_str = '';
		}else{
			$big_seq_str = $big_seq_str.$line;
		}
	}
	push @all_seqs,$big_seq_str;
	shift @all_seqs;
	close IN;
	#print "$#all_seqs +1 seqs obtained ...\n" unless ($input_motif eq 'c');
	
	my $n_pairs = scalar (@all_seqs) /2;
	
	my $found_in_ref = 0;
	my $found_in_tar = 0;
	my $consv_in_tar = 0;
	
	for my $i (0 .. ($n_pairs-1)){
		my $ref_g = $all_seqs[$i*2];
		my $tar_g = $all_seqs[$i*2+1];
		
		#my $found_in_ref = 0;
		#my $consv_in_tag = 0;
		
		my $n1 = length $ref_g;
		my $n2 = length $tar_g;
		
		die("not equal length in ref and tar seq !!!\n$n1 $n2 \n") unless (length $ref_g eq length $tar_g);
		my $L = length $ref_g;
		
		for my $j (0 .. ($L - $l)){
			my $substr_ref = substr($ref_g,$j,$l);
			my $substr_tar = substr($tar_g,$j,$l);
			if(($substr_ref eq $this_motif) && ($substr_tar eq $this_motif)){
				$found_in_ref ++;
				$found_in_tar ++;
				$consv_in_tar ++;
			}elsif ($substr_ref eq $this_motif){
				$found_in_ref ++;
			}elsif ($substr_tar eq $this_motif){
				$found_in_tar ++;	
			}
		}
		
	}


	return ($found_in_ref,$found_in_tar,$consv_in_tar);
}


sub mismatch{
	my $mm1 = shift;
	my $mm2 = shift;
	my @spp1 = split '',$mm1;	
	my @spp2 = split '',$mm2 ;
	my $misn = 0;
	for my $i (0 .. $#spp1){
		$misn ++ unless ($spp1[$i] eq $spp2[$i]);
	}
	return $misn;
}

