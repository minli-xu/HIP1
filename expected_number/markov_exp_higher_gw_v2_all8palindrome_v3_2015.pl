# 1/4/2013 2:41:05 PM modified from analyze_CDS_3RF_baseFreq_dex2012_v3.pl
# 1/18/2013 4:45:09 PM 
# 1/22/2013 4:00:01 PM
# 1/28/2013 3:20:49 PM markov_exp_v1.pl
# 2/13/2013 6:21:05 PM markov higher order
# 2/21/2013 1:50:21 PM v2. correction of rf123, intergenic
# 7/15/2015 4:47:32 PM new genomes

use strict;
use warnings;

my $tag = shift;
my $output = shift;

my $list = 'genome2015_list_5.txt';
my %name;
my %name_NC;
open (LIST,$list) or die('die hard');
while(<LIST>){
	my $line = $_;
	chomp $line;
	my @spline = split "\t",$line;
	$name{$spline[0]} = $spline[1];
	$name_NC{$spline[0]} = $spline[2];
}
close LIST;

open (OUT,">$output") or die('can not crate output file .... ');


my @dirs = <*>;

#my @candidates = ("naz","ama","ava","uca","cyb","cya","cyd","cye","cyf","cyg","mae","ana","npu","sya","sel","syg","syh","syb","syp","syq","syr");
#my @candidates = ("ama","ana","ava","cya","cyb","cyd","cye","cyf","cyg","gvi","mae","naz","npu","pma","pmb","pmc","pmd","pme","pmf","pmg","pmh","pmi","pmj","pmk","pml","sel","sya","syb","syd","sye","syf","syg","syh","syi","sym","syo","syp","syq","syr","uca");
my @one_ple = ('A','C','G','T');
#my @bi_ple;
#my @tri_ple;
#my @tetra_ple;
#my @penta_ple
#my @hexa_ple;
#my @hepta_ple;
#my @octa_ple;
#foreach my $nt1 (@one_ple){
#	foreach my $nt2 (@one_ple){
#		push @bi_ple,$nt1.$nt2;
#		foreach my $nt3 (@one_ple){
#			push @tri_ple,$nt1.$nt2.$nt3;
#			foreach my $nt4 (@one_ple){
#				push @tetra_ple,$nt1.$nt2.$nt3.$nt4;
#				foreach my $nt5 (@one_ple){
#					push @penta_ple,$nt1.$nt2.$nt3.$nt4.$nt5;
#				}
#			}
#		}
#	}
#}

my @motifs;
for my $i1  ($tag .. $tag){
	for my $i2  (0 ..3){
		for my $i3  (0 ..3){
			for my $i4  (0 ..3){
				my $palindrome_motif = $one_ple[$i1].$one_ple[$i2].$one_ple[$i3].$one_ple[$i4].$one_ple[3-$i4].$one_ple[3-$i3].$one_ple[3-$i2].$one_ple[3-$i1];
				push @motifs ,$palindrome_motif;
			}
		}
	}	
}


foreach my $motif (@motifs){

	foreach my $g (@dirs)
	{
		my %hash8mer;
		
		next unless (-d $g);	
		
		my $g_abbr = $name{$g};
		
		#my %params = map { $_ => 1 } @candidates;
		unless(exists($name{$g})) { next; }
	
		chdir ($g);
		
		my @fnas = <*.fna>;
		die ("No fna files found for $g \n") if ((scalar @fnas) == 0);
		
		for my $fna (@fnas){
			my $NC_id = $fna;
			$NC_id =~ s/\.fna//g;
			#print "$NC_id $name_NC{$g}\n";
			next unless ($NC_id eq $name_NC{$g});
			
			open (FNA ,$fna) or die ('can not open fna');
			my %freq;
			my %freq2;
			my %freq3;
			my %freq4;
			my %freq5;
			my %freq6;
			my %freq7;
	
			my $big_str = '';
			while(<FNA>){
				my $line =$_;
				next if ($line =~ /^>/);
				chomp $line;
				$big_str = $big_str.$line;			
			}
			close FNA;
			
			#read the ptt data
			my $ptt = $NC_id.'.ptt';
			open (PTT,$ptt) or die("Can not open the ptt file: $ptt \n");
			my $useless = <PTT>;
			my $genome_size_ptt;
			if ($useless =~ /1\.\.([0-9]+)/)    {
				$genome_size_ptt = $1;
			}else{
				die("Can not get genome size from PTT file : $useless \n");
			}
			$useless = <PTT>;$useless = <PTT>;
			
			my @big_ind = (0) x $genome_size_ptt; ### 0 is intergenic, 1:RF1, 2:RF2, 3:RF3;
			
			while(<PTT>){
				my $line = $_;
				chomp $line;
				my @sp = split "\t",$line;
				my $cord = $sp[0];
				my @c = split '\.\.',$cord;
				next if ($c[0]>$c[1]);
				my $aa_len = $sp[2];
				my @rf1_ind = map -1+0+$c[0]+3*$_,0 .. $aa_len;
				my @rf2_ind = map -1+1+$c[0]+3*$_,0 .. $aa_len;
				my @rf3_ind = map -1+2+$c[0]+3*$_,0 .. $aa_len;
				
				@big_ind[@rf1_ind] = (1) x ($aa_len+1);
				@big_ind[@rf2_ind] = (2) x ($aa_len+1);
				@big_ind[@rf3_ind] = (3) x ($aa_len+1);			
			}
			close PTT;
			
			#end read the ptt data

			
	
	#		my $motif = 'GCGATCGC';
	#		if (($g_abbr eq 'syh') || ($g_abbr eq 'syg')){
	#			$motif = 'GCGATCGC';
	#		}
			my @m = split '',$motif;
			
			#######################
			my @sp_str = split '',$big_str;
			my @region_length = (0,0,0,0);
			for my $i (0 ..$#sp_str){
				last if $i >= $#sp_str-6;
				my $d = $big_ind[$i];
				my $bp2 = $sp_str[$i].$sp_str[$i+1];
				my $bp3 = $sp_str[$i].$sp_str[$i+1].$sp_str[$i+2];
				$freq{$d}{$sp_str[$i]} ++ if defined $freq{$d}{$sp_str[$i]};$freq{$d}{$sp_str[$i]} =1 unless defined $freq{$d}{$sp_str[$i]};
				$freq2{$d}{$bp2} ++ if defined $freq2{$d}{$bp2};$freq2{$d}{$bp2} =1 unless defined $freq2{$d}{$bp2};
				$freq3{$d}{$bp3} ++ if defined $freq3{$d}{$bp3};$freq3{$d}{$bp3} =1 unless defined $freq3{$d}{$bp3};
				$region_length[$d]++;
			}
			my $motif_count = 0;
			for my $i (0 ..$#sp_str){
				last if $i >= $#sp_str-7;
				my $bp8 = $sp_str[$i].$sp_str[$i+1].$sp_str[$i+2].$sp_str[$i+3].$sp_str[$i+4].$sp_str[$i+5].$sp_str[$i+6].$sp_str[$i+7];
				$motif_count ++ if ($bp8 eq $motif);
	
			}		
			my $genome_length = $#sp_str + 1;
			die ("$genome_size_ptt != $genome_length \n") unless ($genome_size_ptt == $genome_length);	
			
			#################
	
			# computing
			
			my @nexp = (0,0,0,0);
			my @nexp2 = (0,0,0,0);
			my @nexp3 = (0,0,0,0);


			for my $d (0 .. 3){
				
				$nexp[$d]  =  $freq{$d}{$m[0]}*$freq{$d}{$m[1]}*$freq{$d}{$m[2]}*$freq{$d}{$m[3]}*$freq{$d}{$m[4]}*$freq{$d}{$m[5]}*$freq{$d}{$m[6]}*$freq{$d}{$m[7]}/($region_length[$d] ** 7);
				
				if (($freq{$d}{$m[1]}*$freq{$d}{$m[2]}*$freq{$d}{$m[3]}*$freq{$d}{$m[4]}*$freq{$d}{$m[5]}*$freq{$d}{$m[6]}) == 0){
					$nexp2[$d] = 0;
				}else{
					$nexp2[$d] =  $freq2{$d}{$m[0].$m[1]}*$freq2{$d}{$m[1].$m[2]}*$freq2{$d}{$m[2].$m[3]}*$freq2{$d}{$m[3].$m[4]}*$freq2{$d}{$m[4].$m[5]}*$freq2{$d}{$m[5].$m[6]}*$freq2{$d}{$m[6].$m[7]} /  ($freq{$d}{$m[1]}*$freq{$d}{$m[2]}*$freq{$d}{$m[3]}*$freq{$d}{$m[4]}*$freq{$d}{$m[5]}*$freq{$d}{$m[6]});
				}
				
				if (($freq2{$d}{$m[1].$m[2]}*$freq2{$d}{$m[2].$m[3]}*$freq2{$d}{$m[3].$m[4]}*$freq2{$d}{$m[4].$m[5]}*$freq2{$d}{$m[5].$m[6]}) == 0){
					$nexp3[$d] = 0;	
				}else{
					$nexp3[$d] =  $freq3{$d}{$m[0].$m[1].$m[2]}*$freq3{$d}{$m[1].$m[2].$m[3]}*$freq3{$d}{$m[2].$m[3].$m[4]}*$freq3{$d}{$m[3].$m[4].$m[5]}*$freq3{$d}{$m[4].$m[5].$m[6]}*$freq3{$d}{$m[5].$m[6].$m[7]} /($freq2{$d}{$m[1].$m[2]}*$freq2{$d}{$m[2].$m[3]}*$freq2{$d}{$m[3].$m[4]}*$freq2{$d}{$m[4].$m[5]}*$freq2{$d}{$m[5].$m[6]});
				}
			}
		
			my $e  = $nexp[0]+$nexp[1]+$nexp[2]+$nexp[3];
			my $e2 = $nexp2[0]+$nexp2[1]+$nexp2[2]+$nexp2[3];
			my $e3 = $nexp3[0]+$nexp3[1]+$nexp3[2]+$nexp3[3];
			
			my $oe  = 'NA';
			my $oe2 = 'NA';
			my $oe3 = 'NA';
			
			$oe  = $motif_count / $e   unless $e == 0;
			$oe2 = $motif_count / $e2  unless $e2 == 0;
			$oe3 = $motif_count / $e3  unless $e3 == 0;

			
			print OUT "$motif\t$g\t$NC_id\t$genome_length\t$motif_count\t$e\t$e2\t$e3\t";
			print OUT "$oe\t$oe2\t$oe3\n";		
			
		}
		
	
		chdir ('..');
	}
	print OUT "\n";

}

close OUT;

# from:
# http://stackoverflow.com/questions/2114185/how-can-i-count-overlapping-substrings-in-perl
sub countnmstr {
    my ($haystack, $needle) = @_;
    my ($first, $rest) = $needle =~ /^(.)(.*)$/;
    return scalar (() = $haystack =~ /(\Q$first\E(?=\Q$rest\E))/g);
}