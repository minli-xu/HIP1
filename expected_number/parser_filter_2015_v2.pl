use strict;
use warnings;

my $oe_cut = 0;
my $n_cut = 0;
my $chi2_cut = 0;
my $f_cut = 0.1;

my @inputs = <markov_exp_*.txt>;

my %nameMap = (
	'Acaryochloris_marina_MBIC11017_uid58167' => 'amr',
	'Anabaena_cylindrica_PCC_7122_uid183339' => 'acy',
	'Anabaena_90_uid179383' => 'anb',
	'Anabaena_variabilis_ATCC_29413_uid58043' => 'ava',
	'Arthrospira_platensis_NIES_39_uid197171' => 'arp',
	'Calothrix_PCC_6303_uid183109' => 'calt',
	'Calothrix_PCC_7507_uid182930' => 'calo',
	'cyanobacterium_UCYN_A_uid43697' => 'cyu',
	'Chroococcidiopsis_thermalis_PCC_7203_uid183002' => 'cthe',
	'Crinalium_epipsammum_PCC_9333_uid183113' => 'cep',
	'Cyanobacterium_aponinum_PCC_10605_uid183340' => 'can',
	'Cyanobacterium_stanieri_PCC_7202_uid183337' => 'csn',
	'Cyanobium_gracile_PCC_6307_uid182931' => 'cgc',
	'Cyanothece_ATCC_51142_uid59013' => 'cyt',
	'Cyanothece_PCC_7424_uid59025' => 'cyc',
	'Cyanothece_PCC_7425_uid59435' => 'cyn',
	'Cyanothece_PCC_7822_uid52547' => 'cyj',
	'Cyanothece_PCC_8801_uid59027' => 'cyp',
	'Cyanothece_PCC_8802_uid59143' => 'cyh',
	'Dactylococcopsis_salina_PCC_8305_uid183341' => 'dsl',
	'Geitlerinema_PCC_7407_uid183007' => 'gei',
	'Gloeobacter_JS_uid225602' => 'glj',
	'Gloeobacter_violaceus_PCC_7421_uid58011' => 'gvi',
	'Gloeocapsa_PCC_7428_uid183112' => 'glp',
	'Halothece_PCC_7418_uid183338' => 'hao',
	'Leptolyngbya_PCC_7376_uid182928' => 'lep',
	'Microcoleus_PCC_7113_uid183114' => 'mic',
	'Microcystis_aeruginosa_NIES_843_uid59101' => 'mar',
	'_Nostoc_azollae__0708_uid49725' => 'naz',
	'Nostoc_punctiforme_PCC_73102_uid57767' => 'npu',
	'Nostoc_PCC_7107_uid182932' => 'nos',
	'Nostoc_PCC_7120_uid57803' => 'ana',
	'Nostoc_PCC_7524_uid182933' => 'nop',
	'Oscillatoria_acuminata_PCC_6304_uid183003' => 'oac',
	'Oscillatoria_nigro_viridis_PCC_7112_uid183110' => 'oni',
	'Pleurocapsa_PCC_7327_uid183006' => 'plp',
	'Prochlorococcus_marinus_AS9601_uid58307' => 'pmb',
	'Prochlorococcus_marinus_MIT_9215_uid58819' => 'pmh',
	'Prochlorococcus_marinus_MIT_9301_uid58437' => 'pmg',
	'Prochlorococcus_marinus_MIT_9303_uid58305' => 'pmf',
	'Prochlorococcus_marinus_MIT_9312_uid58357' => 'pmi',
	'Prochlorococcus_marinus_MIT_9313_uid57773' => 'pmt',
	'Prochlorococcus_marinus_MIT_9515_uid58313' => 'pmc',
	'Prochlorococcus_marinus_NATL1A_uid58423' => 'pme',
	'Prochlorococcus_marinus_NATL2A_uid58359' => 'pmn',
	'Prochlorococcus_marinus_CCMP1375_uid57995' => 'pma',
	'Prochlorococcus_marinus_pastoris_CCMP1986_uid57761' => 'pmm',
	'Pseudanabaena_PCC_7367_uid183004' => 'pseu',
	'Rivularia_PCC_7116_uid182929' => 'riv',
	'Stanieria_cyanosphaera_PCC_7437_uid183115' => 'scs',
	'Synechococcus_elongatus_PCC_6301_uid58235' => 'syc',
	'Synechococcus_elongatus_PCC_7942_uid58045' => 'syf',
	'Synechococcus_CC9311_uid58123' => 'syg',
	'Synechococcus_CC9605_uid58319' => 'syd',
	'Synechococcus_CC9902_uid58323' => 'sye',
	'Synechococcus_JA_2_3B_a_2_13__uid58537' => 'cyb',
	'Synechococcus_JA_3_3Ab_uid58535' => 'cya',
	'Synechococcus_PCC_6312_uid182934' => 'syne',
	'Synechococcus_PCC_7002_uid59137' => 'syp',
	'Synechococcus_PCC_7502_uid183008' => 'synp',
	'Synechococcus_RCC307_uid61609' => 'syr',
	'Synechococcus_WH_7803_uid61607' => 'syx',
	'Synechocystis_PCC_6803_uid189748' => 'syn',
	'Synechocystis_PCC_6803_uid159873' => 'syy',
	'Synechocystis_PCC_6803_uid57659' => 'syz',
	'Synechocystis_PCC_6803_substr__GT_I_uid157913' => 'syt',
	'Synechocystis_PCC_6803_substr__PCC_N_uid159835' => 'syq',
	'Synechocystis_PCC_6803_substr__GT_I_uid158059' => 'sys',
	'Thermosynechococcus_elongatus_BP_1_uid57907' => 'tel',
	'Thermosynechococcus_NK55_uid231517' => 'thn',
	'Trichodesmium_erythraeum_IMS101_uid57925' => 'ter'
);

my %big_oe3;
my %big_chi2;
my %big_f;
my %big_n;

foreach my $f (@inputs){
	open (IN,$f) or die("can not open input file \n")	;
	while(<IN>){
		my $line = $_;
		chomp $line;
		next if $line eq '';
		my @sp = split "\t",$line;
		my $oe3 = $sp[10];
		my $nHIP1 = $sp[4];
		my $eHIP1 = $sp[7];
		my $fHIP1 = $sp[4]/ $sp[3] * 1000;
		my $spp = $nameMap{$sp[1]};
		#die("$spp\n");
		$big_oe3{$sp[0]}{$spp} = $oe3;
		$big_chi2{$sp[0]}{$spp} = ($nHIP1 - $eHIP1)*($nHIP1 - $eHIP1) / $eHIP1;
		$big_f{$sp[0]}{$spp} = $fHIP1;
		$big_n{$sp[0]}{$spp} = $nHIP1;
	}
	close IN;
}

my @all_8mer = sort keys %big_f;
my @all_sp = sort keys (%{$big_f{'AAAATTTT'}});

#$, = "\t";
#print "\t";
#print @all_sp;
#print "\n";

print "#sp\tmotif\tnum\to/e\tchi2\n";
foreach my $m (@all_8mer){
	#print "$m";
	foreach my $sp (@all_sp){
		#print "\t$big_n{$m}{$sp}";
		
		if (($big_chi2{$m}{$sp} >= $chi2_cut) &&  ($big_oe3{$m}{$sp} >= $oe_cut) && ($big_n{$m}{$sp} >= $n_cut) && ($big_f{$m}{$sp} >= $f_cut)){
			print "$sp\t$m\t$big_n{$m}{$sp}\t$big_oe3{$m}{$sp}\t$big_chi2{$m}{$sp}\n";		
		}
	}
	#print "\n"
}