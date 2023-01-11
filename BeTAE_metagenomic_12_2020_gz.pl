#!/usr/bin/perl -w 
#Programme principal
use strict;
use warnings;
use POSIX qw(strftime);
#use Memory::Stats;
# VERSION MG-RAST : ftp://ftp.metagenomics.anl.gov/data/MD5nr/current/
#################################
my $path = $ARGV[0];
my @ins = <$path/*.gz>;
#my @ins = <$path/*.fastq>;
(my $nbsfics) = @ins / 2;
print "\n\n### $nbsfics fics to analyze ###\n\n";
my $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
printf("\n\n### date and time - $datestring ###\n\n");
####### DB storage #######
print "\n\n### Storing Dbs into memory ###\n\n";
#### LCA map ####
my %data_lca;
my $inannots = '/home/terratyv/scratch/Metatranscripto/DBs/md5_lca_map';
print "\n\n\t### Storing DB1 in memory ###\n";
open(DAT,$inannots);
while(<DAT>)
{
        chomp(my $line = $_);
        (my @tab) = split("\t",$line);
        (my $md5) = shift(@tab);
        pop(@tab);
        my $annot = join("\t",@tab);
        $data_lca{$md5}  = $annot;
}
close DAT;
#### Absolute map ####
print "\t### Storing DB2 in memory ###\n";
my %data_tax;
$inannots = '/home/terratyv/scratch/Metatranscripto/DBs/parsed_md5_protein_map';
open(DAT,$inannots);
while(<DAT>)
{
        chomp(my $line = $_);
        (my @tab) = split("\t",$line);
        (my $md5) = shift(@tab);
        my $taxons = join("\t",@tab);
        $data_tax{$md5}  = $taxons;
}
close DAT;
print "\t### Storing DB3 in memory\n\n";
my %tax_line;
my $intax = '/home/terratyv/scratch/Metatranscripto/DBs/parsed.m5nr_v10.taxonomy';
open(DAT,$intax);
while(<DAT>)
{
        chomp(my $line = $_);
        (my @tab) = split("\t",$line);
        (my $tax) = shift(@tab);
        my $taxons = join("\t",@tab);
        $tax_line{$tax}  = $taxons;
}
close DAT;        
####### Main program #######
foreach my $f0_1_gz(@ins)
{
	if( $f0_1_gz =~ /_R1/ )
	{
		my $begin_time = time();
	        (my $sasa) = ($f0_1_gz =~ /([\w]+)_R1/);
		print "PREFIX is $sasa\n";
		my $prepre = 'BeTAE_'.$sasa;
	        if(not -d $prepre)
        	{
			print "\n\n### Generating output directories for $sasa ###\n\n";
		        make_directories($sasa);
			my $pre = 'BeTAE_'.$sasa;
			(my $ref) = $sasa;
			(my $ref1) = $sasa.'_R1';
			(my $ref2) = $sasa.'_R2';
			my $f0_2_gz = $f0_1_gz;
			$f0_2_gz =~ s/_R1/_R2/;
			#(my $f0_1) = ($f0_1_gz =~ /(.*)\.gz/);
			#(my $f0_2) = ($f0_2_gz =~ /(.*)\.gz/);
			(my $ss1) = ($f0_1_gz =~ /\/([^\/]+\.fastq)/);
			(my $ss2) = ($f0_2_gz =~ /\/([^\/]+\.fastq)/);
			#print "\n\n### Gunzip step ###\n\n";
			#system "gunzip $f0_1_gz $f0_2_gz";
			my $f1_1 = $pre.'/SolexaQC/'.$ss1.'.trimmed';
			my $f1_2 = $pre.'/SolexaQC/'.$ss2.'.trimmed';
	                my $SolexaQC = $pre.'/SolexaQC';
		        print "\n\n### Quality trimming ###\n\n";
			solexaqa($f0_1_gz,$f0_2_gz,$SolexaQC);
		        ####
		        #system "rm -rf $f0_1 $f0_2 $f0_1_gz $f0_2_gz";
			#print "\n\n### Gunzip step ###\n\n";
			my $cmdzip = $SolexaQC.'/*.gz';
			system "gunzip $cmdzip";
			print "\n\nremoving short reads (< 75 nt)\n\n";
			length_cutoff1($f1_1,$f1_2);
			my $trimmed1 = $f1_1.'.l75';
			my $trimmed2 = $f1_2.'.l75';	
			#####
			#### KRAKEN rRNA search ####
                	print "\n\n### KRAKEN rRNA search ###\n\n";
 	                my $kraken_report =  $pre.'/rRNA/'.$pre.'.kraken.report';
        	        my $kraken_class = $pre.'/rRNA/'.$pre.'.kraken.classification';
			my $classified = $pre.'/rRNA/'.$pre.'_rRNA';
			my $unclassified = $pre.'/rRNA/'.$pre.'_non_rRNA';
			system "/home/terratyv/scratch/Metatranscripto/kraken2-2.0.8-beta/kraken2 --db /home/terratyv/scratch/Metatranscripto/rRNA_profiles --threads 40 --use-mpa-style --output $kraken_class --report $kraken_report --paired --unclassified-out $unclassified#.fastq --classified-out $classified#.fastq $trimmed1 $trimmed2";
                	clean_kraken($kraken_report, $pre);
			system "rm -rf temp";
			print "\n\n### Fasta conversion of trimmed reads ###\n\n";
			my $out_rep = $pre.'/Assembly_IDBA/'.$pre.'.fasta';
			my $derep1 = $trimmed1.'.derep';
		        my $derep2 = $trimmed2.'.derep';
	        	print "\n\n### Dereplication using the k-mers approach ###\n\n";	        
			derep($trimmed1,$trimmed2);
			system "cat $derep1 $derep2 > $out_rep";
			system "sed -i 's/ /#/g' $out_rep";
                	##### Gene prediction on short reads #####
	                print "\n\n### Gene prediction on short reads ###\n\n";
        	        my $out_FragGeneScan = $pre.'/Gene_Prediction_FragGeneScan/'.$pre;
                	system "/home/terratyv/scratch/Metatranscripto/FragGeneScanPlus/FGS+ -s $out_rep -m 100000 -o $out_FragGeneScan -w 0 -t illumina_10 -p 40";
			###### Fasta cleaning of the FragGenScan output #####
			print "\n\n### Fasta cleaning of the FragGenScan output ###\n\n";
		        my $in_clean = $pre.'/Gene_Prediction_FragGeneScan/'.$pre.'.faa';
		        my $out_clean = $pre.'/Gene_Prediction_FragGeneScan/'.$pre.'_cleaned.faa';
		        clean_that_fas($in_clean,$out_clean);
		        system "rm -rf $in_clean";
		        ###### cd-hit clustering of predicted proteins ######
		        print "\n\n### cd-hit clustering of predicted proteins ###\n\n";
		        my $out_clustering_faa = $pre.'/Gene_Prediction_FragGeneScan/'.$pre.'.C90';
		        my $tempi = 'temp_'.$pre;
			system "/home/terratyv/scratch/Metagenomic/cdhit-4.8.1/cd-hit -i $out_clean -o $out_clustering_faa -M 0 -T 0 -c 0.9 -d 200 > $tempi";
		        print "\n\n### Cluster abundance calculation for protein predictions ###\n\n";
		        my $in_abundance_faa = $pre.'/Gene_Prediction_FragGeneScan/'.$pre.'.C90.clstr';
		        (my %abd_faa) = abundance($in_abundance_faa);
		        ###### DIAMOND on MD5nr ######
		        print "\n\n### Diamond search on MD5nr ###\n\n";
		        my $outfile = $pre.'/Similarity_search/Diamond_against_md5nr_'.$pre;
			system "/home/terratyv/scratch/Metatranscripto/diamond/diamond blastp -d /home/terratyv/scratch/Metatranscripto/M5nr/md5nr.dmnd -q $out_clustering_faa -o $outfile -e 1E-05 --threads 40 > $tempi";
			system "rm -rf $tempi";
			###### Getting hits for MD5nr ######
		        print "\n\n### Best hits retrieval for MD5nr ###\n\n";
		        my $evalue = 1E-05;
		        my ($res,$correction) = best_hits($outfile,$evalue);
		        #############################
		        ##### Taxonomic profiles ####
		        #############################
		        print "\n\n### Taxonomic profile of MD5nr hits using LCA ###\n\n";
		        my $begin = time();
		        my ($nhits,$nannot) = lca_request_hash($res,\%abd_faa,$correction,$pre,$ref,\%data_lca);
		        my $end = time() - $begin;
		        print "\n\n#### LCA profiles ####\n\nNb hits : $nhits Nb annots : $nannot\nUsed $end for search ###\n\n";
		        ##### Taxonomic profile of MD5nr hits from absolute taxonomy using SQlite #####
		        print "\n\n### Taxonomic profile of MD5nr hits using Absolute names ###\n\n";
		        $begin = time();
	 		($nhits,$nannot) = absolute_request_hash($res,\%abd_faa,$correction,$pre,$ref,\%data_tax,\%tax_line);
	 		print "\n\n#### Absolute profiles ####\n\nNb hits : $nhits Nb annots : $nannot\nUsed $end for SQlite search ###\n\n";
			#############################
		        #############################
		        ##### Functional profiles ####
		        #############################
		        #print "### Generating Subsystems profiles using SQLite ###\n\n";
		        my $table = 'SUBSYSTEMS';
		        my $nb_fields = 4;
		        ######
		        print "\n\n### Generating Subsystems profiles###\n\n";
		        $begin = time();
		        ($nhits,$nannot) = ontology_hash($res,\%abd_faa,$correction,$table,$nb_fields,$pre,$ref);
			print "\n\n#### Subsystems profiles ####\n\nNb hits : $nhits Nb annots : $nannot\nUsed $end for HASH search ###\n\n";
			######				
		        print "\n\n### Generating COGs profiles ###\n\n";
		        $table = 'COGs';
		        $nb_fields = 3;
		        ($nhits,$nannot) = ontology_hash($res,\%abd_faa,$correction,$table,$nb_fields,$pre,$ref);
			print "\n\n#### COGs profiles ####\n\nNb hits : $nhits Nb annots : $nannot\nUsed $end for SQlite search ###\n\n";
		        ##### KOs profiles using SQlite #####
		        print "### Generating KOs profiles ###\n\n";
		        $table = 'KOs';
		        $nb_fields = 4;
		        ($nhits,$nannot) = ontology_hash($res,\%abd_faa,$correction,$table,$nb_fields,$pre,$ref);
			print "\n\n#### KOs profiles ####\n\nNb hits : $nhits Nb annots : $nannot\nUsed $end for SQlite search ###\n\n";
			my $end_time = time();
			my $time_used =  sprintf("%.1f", ($end_time - $begin_time) / 60);
			print "##################################\n\nUsed $time_used minutes\n\n##################################\n\n";
			#### Simple Names profile ####
			#$table = 'SIMPLE NAMES';
		        #$nb_fields = 1;
		        #print "### Generating Simple Names ###\n\n";
		        #$begin = time();
		        #($nhits,$nannot) = ontology_hash($res,\%abd_faa,$correction,$table,$nb_fields,$pre,$ref);
			#print "\n\n#### Subsystems profiles ####\n\nNb hits : $nhits Nb annots : $nannot\nUsed $end for HASH search ###\n\n";
			print "\n\n### Cleaning inputs ###\n\n";
			#system "rm -rf $f0_1 $f0_2"; 
			print "\n\n###########################################################\n###########################################################\n\n";
		}
	}
}
exit;
#################################
####### Subroutines       #######
#################################
#################################
sub fastq2fasta
{
	(my @data) = @_;
	open(OUT,">$data[2]");
	for(my $i =0;$i <2 ;$i ++)
	{
		open(IN,$data[$i]);
		my $cl =0;
		while(<IN>)
		{
			$cl ++;
			if($cl ==1)
			{	
				chomp(my $line = $_);
				$line =~ s/^\@/>/;
				print OUT "$line\n";
			}
			elsif($cl ==2)
			{
				print  OUT "$_";
			}
			elsif($cl ==4)
			{
				$cl = 0;
			}
		}
		close IN;
	}		
	close OUT;
}
###############
###############
###############
###############
sub clean_kraken
{
	(my @data) = @_;
	open(IN,$data[0]);
	my $pre = $data[1];
	my $out = $pre.'/rRNA/Matrix_rRNA_'.$pre.'.out';
	open(OUT,">$out");
	my $tot = 0;
	my $res;
	my $unclass = 0;;
	while(<IN>)
	{
		(my @tab) = split("\t",$_);
		#d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Rhizobiales|f__Methylobacteriaceae|s__Methylobacteriaceae bacterium KVD-1982-08	7955
		#d__Viroids|f__Pospiviroidae|g__Hostuviroid|s__Hop
		my $domain = 'Unclassified';
		my $phylum = 'Unclassified';
		my $class = 'Unclassified';
		my $order = 'Unclassified';
		my $family = 'Unclassified';
		my $genus = 'Unclassified';
		my $species = 'Unclassified';
		if($tab[0] =~ /s__/)
		{
			if($tab[0] =~ /d__([^\|]+)\|/)
			{
				$domain = $1;
			}
			if($tab[0] =~ /p__([^\|]+)\|/)
			{
				$phylum = $1;
			}
			if($tab[0] =~ /c__([^\|]+)\|/)
			{
				$class = $1;
			}
			if($tab[0] =~ /o__([^\|]+)\|/)
			{
				$order = $1;
			}
			if($tab[0] =~ /f__([^\|]+)\|/)
			{
				$family = $1;
			}
			if($tab[0] =~ /g__([^\|]+)\|/)
			{
				$genus = $1;
			}
			($species) = ($tab[0] =~ /s__([^\|\t]+)/)
		}
		my $tax = "$domain\t$phylum\t$class\t$order\t$family\t$genus\t$species";
		my $nb = 0;
		while($tax =~ m/Unclassified/ig)
		{
			$nb ++;
		}
		if($nb == 7)
		{
			$unclass += $tab[1];
		}
		else
		{
			print OUT "$tax\t$tab[1]";
		}
	}
	print OUT "Unclassified\tUnclassified\tUnclassified\tUnclassified\tUnclassified\tUnclassified\tUnclassified\t$unclass";
	close OUT;
}
#####################################################
#####################################################
sub get_fasta_from_clustering
{
	(my @data) = @_;
	#get_fasta_from_clustering($outfile_fna,$in_abundance_fna, $fas_rRNA_pre_screen);
	my $name;
	my %fasta;
	open(FAS,$data[0]);
	while(<FAS>)
	{
		if($_ =~ /^>/)
		{
			chomp($name = $_);
		}
		else
		{
			chomp($fasta{$name} .= $_);
		}
	}
	close FAS;
	open(CLSTR,$data[1]);
	open(OUT,">$data[2]");
	while(<CLSTR>)
	{
		if($_ =~ /(>.*)\.{3} \*/)
		{
			print OUT "$1\n$fasta{$1}\n";
		}
	}
	close CLSTR;
	close OUT;
}
#################################
#################################
#################################
sub interleaved_fasta
{
	(my @data) = @_;
	my $in1 = $data[0];
	my $in2 = $data[1];
	my $out = $data[2];
	open(IN1,$in1);
	open(IN2,$in2);
	open(OUT,">$out");
	my $name;
	my %f1;
	my $c1 = 0;
	while(<IN1>)
	{
	        if($_ =~ /^>/)
                {
                        chomp($name = $_);
                        $c1 ++;
                }
                else
                {
                        chomp($f1{$name} .= $_);
                }
     
        }
        close IN1;
        my %f2;
        my $c2 = 0;
	while(<IN2>)
	{
	        if($_ =~ /^>/)
                {
                        chomp($name = $_);
                        $c2 ++;
                }
                else
                {
                        chomp($f2{$name} .= $_);
                }
     
        }
        close IN2;
        if($c1 >= $c2)
        {
                (my @sort) = sort{$a cmp $b} keys (%f1);
                foreach my $e1(@sort)
                {
                        my $e2 = $e1;
                        $e2 =~  s/ 1:/ 2:/;
                        if($f2{$e2})
                        {
                                print OUT "$e1\n$f1{$e1}\n$e2\n$f2{$e2}\n";
                        }
                }
        }
        else
        {
                (my @sort) = sort{$a cmp $b} keys (%f2);
                foreach my $e2(@sort)
                {
                        my $e1 = $e2;
                        $e1 =~  s/ 2:/ 1:/;
                        if($f1{$e1})
                        {
                                print OUT "$e1\n$f1{$e1}\n$e2\n$f2{$e2}\n";
                        }
                }
        }
        %f1 = ();
        %f2 = ();
}
#################################
#################################
sub position_biases
{
	(my @data) = @_;
	my $in = $data[0];
	(my $tit) = ($in =~ /(\d+)_R/);
	my $out = $data[1];
	my $out2 = $data[2];
	my $pre = $data[3];
	open(IN,$in);
	my $pos;
	my $tot = 0;
	my $c = 0;
	while(<IN>)
	{
	        if($c == 0)
                {
                        $c ++;
                }
                elsif($c == 1)
                {
                        chomp(my $seq = $_);
                        $tot = length($seq);
                        last;
                }
        }
        close IN;
        for(my $i = 0; $i < $tot; $i ++)
        {
                $$pos{$i}{'A'} = 0;
                $$pos{$i}{'C'} = 0;
                $$pos{$i}{'T'} = 0;
                $$pos{$i}{'G'} = 0;
                $$pos{$i}{'N'} = 0;
        }
        $c = 0;
        open(IN,$in);
	while(<IN>)
	{
	        if($c == 0)
                {
                        $c ++;
                }
                elsif($c == 1)
                {
                        chomp(my $seq = $_);
                        (my @posi) = split("",$seq);
                        my $count = 0;
                        #print "\n\n$seq\n\n";
                        foreach my $e(@posi)
                        {
                                $$pos{$count}{$e} ++;
                                #print "\t$count : $e\n";
                                $count ++;
                        }
                        $c ++;
                        #sleep(10);
                }
                elsif($c == 2)
                {
                        $c ++;
                }
                elsif($c == 3)
                {
                        $c = 0;
                }
        }
        close IN;
        open(OUT,">$out");
        print OUT "Position\tA\tC\tT\tG\tN\n";
        my $line;
        for(my $i = 0; $i < $tot; $i ++)
        {
                print OUT "$i\t$$pos{$i}{'A'}\t$$pos{$i}{'C'}\t$$pos{$i}{'T'}\t$$pos{$i}{'G'}\t$$pos{$i}{'N'}\n";
                $line .= "$$pos{$i}{'A'}\n$$pos{$i}{'C'}\n$$pos{$i}{'T'}\n$$pos{$i}{'G'}\n$$pos{$i}{'N'}\n";
        }
        close OUT;        
        open(OUT2,">$out2");
        $line =~ s/\n$//;
        print OUT2 "$line\n";
        close OUT2;
        my $rsc = $pre.'/SolexaQC/Rscript.position.biases.R';
        open(R,">$rsc");
        my $outpng =  $pre.'/SolexaQC/'.$tit.'.png';
        my $cline =  'library(ggplot2)'."\n".'library("RColorBrewer")'."\n".'png("'.$outpng.'",width=2048,height=2048)'."\n".'frequency<-read.csv("'.$out2.'", header=F)'."\n".'base=rep(c("A","C","T","G","X"),'.$tot.')'."\n".'position=gl('.$tot.',5)'."\n".'data=data.frame(position,base,frequency)'."\n".'col= c(brewer.pal(n = 4, name = "PiYG"),"black")'."\n".'ggplot(data, aes(fill=base, y=frequency, x=position)) + geom_bar( stat="identity", position="fill") + scale_fill_manual(values=col) + theme(axis.text.x = element_text(face="bold", size = 30,angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(face="bold", size = 30, hjust = 1, vjust = 0.5), axis.title.x = element_text(face="bold", size=40) , axis.title.y = element_text(face="bold",  size=40), plot.title = element_text(size=70, face="bold",hjust = 0.5, vjust = 0.5), legend.title=element_blank(), legend.text = element_text(size=70, face="bold")) + ggtitle("'.$tit.'") + scale_x_discrete(breaks=seq(0,'.$tot.',5))'."\n".'dev.off()';
        print R "$cline\n";
        close R;
        system "Rscript  $rsc";      
}
#################################
#################################
#################################
#################################
sub shannon_diversity_index
{
	(my @data) = @_;
	my $in = $data[0];
	my $out = $data[1];
	open(IN,$in);
	open(OUT,">$out");
	
}
#################################
#################################
sub solexaqa
{
        (my @data) = @_;
	my $in1 = $data[0];
	my $in2 = $data[1];
	my $out = $data[2];
        my @pid;
        my $nbtot = 2;
	my @bioms_thread;
	push(@bioms_thread,$in1,$in2);
	my $c = -1;
	for(my $i=0; $i < 2; $i++)
	{
		$c ++;
		push(@pid,my $pid=fork);
		unless($pid)
		{
			system "/home/terratyv/scratch/Metatranscripto/SolexaQA++/SolexaQA++ dynamictrim $bioms_thread[$c] -d $out";
                       exit(0);
                }
	}
	foreach (@pid) 
       {
	        waitpid($_,0);
        }
}
#################################
#################################
sub derep
{
	(my @data) = @_;
	my $in1 = $data[0];
	my $in2 = $data[1];
	my @pid;
        my $nbtot = 2;
	my @bioms_thread;
	push(@bioms_thread,$in1,$in2);
	my $c = -1;
	for(my $i=0; $i < 2; $i++)
	{
		$c ++;
		push(@pid,my $pid=fork);
		unless($pid)
		{
			my $derep = $bioms_thread[$c].'.derep';
			open(IN,$bioms_thread[$c]);
			open(OUT, ">$derep");
                	my %kmers;
                        my %fas;
                        my %refs;
                        my $name;
                        my $c = 0;
                        my $init = 0;
                        while(<IN>)
	                {
	                        if($c == 0)
                                {
                                        chomp($name = $_);
                                        $c ++;
                                        $init ++;
                                }
                                elsif($c == 1)
                                {
                                        chomp(my $seq = $_);
                                        $fas{$name} = $seq;
                                        (my $pre) = substr($seq,0,19);
                                        $kmers{$pre} ++;
                                        $refs{$pre} .= $name.'--';
                                        $c ++;
                                }
                                elsif($c == 2)
                                {
                                        $c ++;
                                }
                                elsif($c == 3)
                                {
                                        $c = 0;
                                }
                        }
                        close IN;
                        my $select = 0;
                        (my @sort) = sort{$kmers{$b} <=> $kmers{$a}} keys (%kmers);
                        # Optimiser en prenant la plus grande sequence #
                        foreach my $e(@sort)
                        {
                                $select ++;
                                (my @seqs) = split('--',$refs{$e});
                                my $nn = $seqs[0];
                                $nn =~ s/^@/>/;
                                print OUT "$nn\n$fas{$seqs[0]}\n";
                        }
                        my $ratio = sprintf("%.4f", $select/$init) * 100;
                        my $rep = sprintf("%.2f",100 - $ratio);
                        print "\t\t$ratio % of trimmed reads selected for further analysis\n\t\t$rep % of trimmed reads are identified as artificial duplicates\n\n";
                        %kmers = ();
                        %fas = ();
                        %refs = ();
                        exit(0);
                }
	}
	foreach (@pid) 
        {
	        waitpid($_,0);
        }
}	
##################################
##################################
sub ontology_hash
{
	(my @data) = @_;
        my $res = $data[0];
        my $abd = $data[1];
        my $cor = $data[2];
        my $table = $data[3];
        my $nb_fields = $data[4];
        my $pre = $data[5];
        my $ref = $data[6];
        my %data;
        my $inannots = '/home/terratyv/scratch/Metatranscripto/DBs/'.$table.'_complete_DB';
        open(DAT,$inannots);
        while(<DAT>)
        {
                chomp(my $line = $_);
                (my @tab) = split("\t",$line);
                (my $md5) = shift(@tab);
                my $annot = join("\t",@tab);
                $data{$md5}  = $annot;
        }
        close DAT;
	############################################################################################################
	my $profile_ortho;
	my $profile_ortho_corrected;
	my %profile_ortho_comp;
	my %profile_ortho_corrected_comp;
	my $nb_hits = 0;
	my $nb_annot = 0;
	foreach my $k(keys(%$res))
        {
                my $factor = $$abd{$k};
                my $factor_corrected = $$abd{$k} / $$cor{$k};
                $nb_hits += $factor;
                foreach my $md5(keys(%{$res -> {$k}}))
                {
			if($data{$md5})
			{
			        (my @infs) = split("\t",$data{$md5});
			        for(my $i = 0; $i < 2; $i ++)
			        {
			                shift(@infs);
		                }
			        (my $ter) = join("\t",@infs);
			        if($ter ne '')
			        {
			                $nb_annot += $factor;
			                $profile_ortho_comp{$ter} += $factor;
			                $profile_ortho_corrected_comp{$ter} += $factor_corrected;
		                }
	                }
		}
        }
        ############################################################################################################
	my $outcomp = $pre.'/Functions/'.$table.'/matrix_'.$table.'_best_hits_'.$pre.'.out';
	open(OUTCOMP,">$outcomp");
	if($table =~ /KOs|SUBSYSTEMS/)
	{
	        print OUTCOMP "Level1\tLevel2\tLevel3\tLevel4\tcount\n";
        }
        else
        {
              print OUTCOMP "Level1\tLevel2\tLevel3\tcount\n";
        }  
	(my @sorta) = sort{$a cmp $b} keys(%profile_ortho_comp);
	foreach my $e(@sorta)
	{
	       print OUTCOMP "$e\t$profile_ortho_comp{$e}\n";
        }
	close OUTCOMP;
	############################################################################################################
	my $outcompcorr = $pre.'/Functions/'.$table.'/matrix_'.$table.'_best_hits_corrected_'.$pre.'.out';
	open(OUTCOMPCORR,">$outcompcorr");
	if($table =~ /KOs|SUBSYSTEMS/)
	{
	        print OUTCOMPCORR "Level1\tLevel2\tLevel3\tLevel4\tcount\n";
        }
        else
        {
              print OUTCOMPCORR "Level1\tLevel2\tLevel3\tcount\n";
        }  
	(@sorta) = sort{$a cmp $b} keys(%profile_ortho_corrected_comp);
	foreach my $e(@sorta)
	{
	        print OUTCOMPCORR "$e\t$profile_ortho_corrected_comp{$e}\n";
        }
	close OUTCOMPCORR;	
        ############################################################################################################
	my $profile_ortho_representative;
	my %profile_representative_comp;
	foreach my $k(keys(%$res))
        {
                my $factor = $$abd{$k};
                foreach my $md5(keys(%{$res -> {$k}}))
                {
			if($data{$md5})
			{
			        (my @infs) = split("\t",$data{$md5});
			        for(my $i = 0; $i < 2; $i ++)
			        {
			                shift(@infs);
		                }
			        (my $ter) = join("\t",@infs);
			        if($ter ne '')
			        {
			                $profile_representative_comp{$ter} += $factor;
		                }
	                }
	       	        last;
	       	}
        }
	############################################################################################################
	$outcomp = $pre.'/Functions/'.$table.'/matrix_'.$table.'_representative_'.$pre.'.out';
	open(OUTCOMP,">$outcomp");
	if($table =~ /KOs|SUBSYSTEMS/)
	{
	        print OUTCOMP "Level1\tLevel2\tLevel3\tLevel4\tcount\n";
        }
        else
        {
              print OUTCOMP "Level1\tLevel2\tLevel3\tcount\n";
        }  
	(@sorta) = sort{$a cmp $b} keys(%profile_representative_comp);
	foreach my $e(@sorta)
	{
	        print OUTCOMP "$e\t$profile_representative_comp{$e}\n";
        }
	close OUTCOMP;
	############################################################################################################
	%data = ();
	$profile_ortho = ();
	$profile_ortho_corrected = ();
	%profile_ortho_comp = ();
	%profile_ortho_corrected_comp = ();
	return($nb_hits,$nb_annot);
}
####################################################
####################################################
sub lca_request_hash
{
	(my @data) = @_;
        my $res = $data[0];
        my $abd = $data[1];
        my $cor = $data[2];
        my $pre = $data[3];
        my $ref = $data[4];
        my $data_lca = $data[5];
        #my $inannots = '/media/yves/Xtra/MG-RAST_like/DEF/usefull_DBs/ftp.metagenomics.anl.gov/FINAL/md5_lca_map';
        #open(DAT,$inannots);
        #while(<DAT>)
        #{
        #        chomp(my $line = $_);
        #        (my @tab) = split("\t",$line);
        #        (my $md5) = shift(@tab);
        #        pop(@tab);
        #        my $annot = join("\t",@tab);
        #        $data{$md5}  = $annot;
        #}
        #close DAT;
	############################################################################################################
	my $profile_lca;
	my $profile_lca_corrected;
	my %profile_lca_comp;
	my %profile_lca_corrected_comp;
	my $nb_hits = 0;
	my $nb_annot = 0;
	foreach my $k(keys(%$res))
        {
                my $factor = $$abd{$k};
                my $factor_corrected = $$abd{$k} / $$cor{$k};
                $nb_hits += $factor;
                foreach my $md5(keys(%{$res -> {$k}}))
                {
			if($$data_lca{$md5})
			{
			        (my @infs) = split("\t",$$data_lca{$md5});
			        (my $ter) = join("\t",@infs);
			        if($ter ne '')
			        {
			                $nb_annot += $factor;
			                $profile_lca_comp{$ter} += $factor;
			                $profile_lca_corrected_comp{$ter} += $factor_corrected;
		                }
		        }	
	       	}
        }
        ############################################################################################################
        ############################################################################################################
	my $outcomp = $pre.'/Taxonomy/LCA/matrix_best_hits_LCA_'.$pre.'.out';
	open(OUTCOMP,">$outcomp");
	print OUTCOMP "domain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\tcount\n";
	(my @sorta) = sort{$a cmp $b} keys(%profile_lca_comp);
	foreach my $e(@sorta)
	{
	       print OUTCOMP "$e\t$profile_lca_comp{$e}\n";
        }
	close OUTCOMP;
	############################################################################################################
	my $outcompcorr = $pre.'/Taxonomy/LCA/matrix_best_hits_corrected_LCA_'.$pre.'.out';
	open(OUTCOMPCORR,">$outcompcorr");
	print OUTCOMPCORR "domain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\tcount\n";
	(@sorta) = sort{$a cmp $b} keys(%profile_lca_corrected_comp);
	foreach my $e(@sorta)
	{
	        print OUTCOMPCORR "$e\t$profile_lca_corrected_comp{$e}\n";
        }
	close OUTCOMPCORR;	
	############################################################################################################
	my $profile_lca_representative;
	my %profile_lca_representative_comp;
	print "\n\n\n\n";
	foreach my $k(keys(%$res))
        {
                my $factor = $$abd{$k};
                foreach my $md5(keys(%{$res -> {$k}}))
                {
			if($$data_lca{$md5})
			{
			        (my @infs) = split("\t",$$data_lca{$md5});
			        (my $ter) = join("\t",@infs);
			        if($ter ne '')
			        {
			                $profile_lca_representative_comp{$ter} += $factor;
		                }
		        }
	       	        last;
	       	}
        }
	############################################################################################################
	############################################################################################################
	$outcomp = $pre.'/Taxonomy/LCA/matrix_representative_LCA_'.$pre.'.out';
	open(OUTCOMP,">$outcomp");
	print OUTCOMP "domain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\tcount\n";
	(@sorta) = sort{$a cmp $b} keys(%profile_lca_representative_comp);
	foreach my $e(@sorta)
	{
	        print OUTCOMP "$e\t$profile_lca_representative_comp{$e}\n";
        }
	close OUTCOMP;	
	#%data = ();
        $profile_lca = ();
	$profile_lca_corrected = ();
	%profile_lca_comp = ();
	%profile_lca_corrected_comp = ();
	$profile_lca_representative =();
	%profile_lca_representative_comp =();
        return($nb_hits,$nb_annot);
}
##################################
#################################
#################################
sub absolute_request_hash
{
	(my @data) = @_;
        my $res = $data[0];
        my $abd = $data[1];
        my $cor = $data[2];
        my $pre = $data[3];
        my $ref = $data[4];
	my $data_tax = $data[5];
	my $tax_line = $data[6];
        #print "\t\tStoring MD5 map in memory\n";
        #my $inannots = '/media/yves/Xtra/MG-RAST_like/DEF/usefull_DBs/ftp.metagenomics.anl.gov/FINAL/parsed_md5_protein_map';
        #my %data;
        #open(DAT,$inannots);
        #while(<DAT>)
        #{
        #        chomp(my $line = $_);
        #        (my @tab) = split("\t",$line);
        #        (my $md5) = shift(@tab);
        #        my $taxons = join("\t",@tab);
        #        $data{$md5}  = $taxons;
        #}
        #close DAT;
        #print "\t\tStoring Taxonomy map in memory\n";
        #my %tax_line;
        #my $intax = '/media/yves/Xtra/MG-RAST_like/DEF/usefull_DBs/ftp.metagenomics.anl.gov/FINAL/parsed.m5nr_v10.taxonomy';
        #open(DAT,$intax);
        #while(<DAT>)
        #{
        #        chomp(my $line = $_);
        #        (my @tab) = split("\t",$line);
        #        (my $tax) = shift(@tab);
        #        my $taxons = join("\t",@tab);
        #        $tax_line{$tax}  = $taxons;
        #}
        #close DAT;       
	############################################################################################################
	my %profile_tax_comp;
	my %profile_tax_corrected_comp;
	my $nb_hits = 0;
	my $nb_annot = 0;
	foreach my $k(keys(%$res))
        {
                my $factor = $$abd{$k};
                my $factor_corrected = $$abd{$k} / $$cor{$k};
                $nb_hits += $factor;
                foreach my $md5(keys(%{$res -> {$k}}))
                {
			if($$data_tax{$md5})
	                {
		                $nb_annot ++;
	                        (my @taxons) = split(" ",$$data_tax{$md5});
				my $nb_taxons = @taxons;
	                        $factor_corrected = $factor_corrected  / $nb_taxons;
	                        foreach my $t(@taxons)
	        		{
	        			$profile_tax_comp{$$tax_line{$t}} += $factor;
				        $profile_tax_corrected_comp{$$tax_line{$t}} += $factor_corrected;
			        }
		        }
	        }
        }
        ############################################################################################################
        ############################################################################################################
	my $outcomp = $pre.'/Taxonomy/Absolute/matrix_best_hits_Absolute_Taxonomy_M5nr'.$pre.'.out';
	open(OUTCOMP,">$outcomp");
	print OUTCOMP "domain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\tcount\n";
	(my @sorta) = sort{$a cmp $b} keys(%profile_tax_comp);
	foreach my $e(@sorta)
	{
	       print OUTCOMP "$e\t$profile_tax_comp{$e}\n";
        }
	close OUTCOMP;
	############################################################################################################
	my $outcompcorr = $pre.'/Taxonomy/Absolute/matrix_best_hits_corrected_Absolute_Taxonomy_M5nr'.$pre.'.out';
	open(OUTCOMPCORR,">$outcompcorr");
	print OUTCOMPCORR "domain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\tcount\n";
	(@sorta) = sort{$a cmp $b} keys(%profile_tax_corrected_comp);
	foreach my $e(@sorta)
	{
	        print OUTCOMPCORR "$e\t$profile_tax_corrected_comp{$e}\n";
        }
	close OUTCOMPCORR;	
	############################################################################################################
	my %profile_tax_representative_comp;
	foreach my $k(keys(%$res))
        {
                my $factor = $$abd{$k};
                foreach my $md5(keys(%{$res -> {$k}}))
                {
			if($$data_tax{$md5})
	                {
		                (my @taxons) = split(" ",$$data_tax{$md5});
				my $nb_taxons = @taxons;
	                        $factor = $factor  / $nb_taxons;
	                        foreach my $t(@taxons)
	        		{
	        			$profile_tax_representative_comp{$$tax_line{$t}} += $factor;
			        }
		        }
		        last;
	        }
        }
	############################################################################################################
	############################################################################################################
	$outcompcorr = $pre.'/Taxonomy/Absolute/matrix_representative_hit_Absolute_Taxonomy_M5nr_'.$pre.'.out';
	open(OUTCOMPCORR,">$outcompcorr");
	print OUTCOMPCORR "domain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\tcount\n";
	(@sorta) = sort{$a cmp $b} keys(%profile_tax_representative_comp);
	foreach my $e(@sorta)
	{
	        print OUTCOMPCORR "$e\t$profile_tax_representative_comp{$e}\n";
        }
	close OUTCOMPCORR;	
	############################################################################################################
	%profile_tax_comp =();
	%profile_tax_corrected_comp =();
	%profile_tax_representative_comp= ();
	#%data = ();
	#%tax_line = ();
        return($nb_hits,$nb_annot);
}
##################################
##################################
##################################
##################################
sub get_rRNA_hits
{
	(my @data) = @_;
	my $in = $data[0];
	my $eval = '1E-05';
	my %res;
	foreach my $inbl(@$in)
	{
		open(IN,$inbl);
		#print "\t\t\tParsing BLAT $inbl output\n";
		while(<IN>)
		{
			(my @tab) = split("\t",$_);
			if( $tab[10] <= $eval )
			{
				$res{$tab[0]} ++;
			}
		}
		close IN;
	}
	my $name;
	my %fasta;
	print "\t\t\tStoring FAS in memory\n";
	open(FAS,$data[1]);
	while(<FAS>)
	{
	        chomp(my $line = $_);
	        if($line =~ /^>(.*)/)
	        {
	                $name = $1;
	                $name =~ s/\s+/#/;
                }
                else
                {
                        $fasta{$name} .= $line;
                }
        }
        close FAS;
        my $outrRNA = $data[2];
        print "\t\t\tGenerating file with pre-screen rRNA\n\n";
        open(OUTRNA,">$outrRNA");
        foreach my $k(keys(%res))
        {
                #my $rename = $k;
                #$rename =~ s/#/ /;
                print OUTRNA ">$k\n$fasta{$k}\n";
        }
        close OUTRNA;
	%res=();
        %fasta=();
}
##################################
##################################
##################################
##################################
sub get_rRNA_hits2
{
	(my @data) = @_;
	my $in = $data[0];
	my $eval = '1E-05';
	my %res;
	open(IN,$in);
	while(<IN>)
	{
		(my @tab) = split("\t",$_);
		if( $tab[10] <= $eval )
		{
			$res{$tab[0]} ++;
		}
	}
	close IN;
	my $name;
	my %fasta;
	print "\t\t\tStoring FAS in memory\n";
	open(FAS,$data[1]);
	while(<FAS>)
	{
	        chomp(my $line = $_);
	        if($line =~ /^>(.*)/)
	        {
	                $name = $1;
	                $name =~ s/\s+/#/;
                }
                else
                {
                        $fasta{$name} .= $line;
                }
        }
        close FAS;
        my $outrRNA = $data[2];
        print "\t\t\tGenerating file with pre-screen rRNA\n\n";
        open(OUTRNA,">$outrRNA");
        foreach my $k(keys(%res))
        {
                #my $rename = $k;
                #$rename =~ s/#/ /;
                print OUTRNA ">$k\n$fasta{$k}\n";
        }
        close OUTRNA;
	%res=();
        %fasta=();
}
##################################
##################################
sub best_hits
{
	(my @data) = @_;
	my $in = $data[0];
	my $eval = $data[1];
	open(IN,$in);
	my $temp;
	my $res;
	my %max;
	while(<IN>)
	{
		(my @tab) = split("\t",$_);
		if($tab[10] <= $eval )
		{
			$$temp{$tab[0]}{$tab[10]}{$tab[1]} ++;
		}
	}
	close IN;
	foreach my $k(keys(%$temp))
	{
	        (my @soreval) = sort{$$temp{$k}{$b} <=> $$temp{$k}{$a}} keys(%{$temp ->{$k}});
	        foreach my $t(keys(%{$temp->{$k}->{$soreval[0]}}))
	        {
	               $$res{$k}{$t} ++;
	               $max{$k} ++;
	        }
        }
        $temp = undef;
        return($res,\%max);
}	


#################################
#################################
sub matches
{
	(my @data) = @_;
	my $in = $data[0];
	my $eval = $data[1];
	open(IN,$in);
	my $res;
	while(<IN>)
	{
		(my @tab) = split("\t",$_);
		if( $tab[10] <= $eval )
		{
			$$res{$tab[0]}{$tab[1]} ++;
		}
	}
	close IN;
	return($res);
}
#################################
#################################
sub length_cutoff1
{
        (my @data) = @_;
        my $in1 = $data[0];
	my $in2 = $data[1];
	my $out1 = $in1.'.l75';
	my $out2 = $in2.'.l75';
	my @pid;
        open(IN1,$in1);
        open(OUT1,">$out1");
        my $cl =0;
        my @block;
        print "\t\t$in1 lenght cleaning\n\n";
	while(<IN1>)
        {
                $cl ++;
                if($cl < 4)
                {
                        push(@block,$_);
                }
                else
                {
                        push(@block,$_);
                        $cl = 0;
                        if(length($block[1]) >= 75)
                        {
                                print OUT1 "@block";
                        }
                        @block = ();
                }
        }
        close IN1;
        close OUT1;
        open(IN2,$in2);
        open(OUT2,">$out2");
        $cl = 0;
	print "\t\t$in2 lenght cleaning\n\n";
        while(<IN2>)
        {
                $cl ++;
                if($cl < 4)
                {
                        push(@block,$_);
                }
                else
                {
                        push(@block,$_);
                        $cl = 0;
                        if(length($block[1]) >= 75)
                        {
                                print OUT2 "@block";
                        }
                        @block = ();
                }
        }
        close IN2;
        close OUT2;
}
#################################
#################################
#################################
sub length_cutoff
{
        (my @data) = @_;
        my $in1 = $data[0];
	my $in2 = $data[1];
        my @pid;
        my $nbtot = 2;
	my @bioms_thread;
	push(@bioms_thread,$in1,$in2);
	my $c = -1;
	for(my $i=0; $i < 2; $i++)
	{
		$c ++;
		push(@pid,my $pid=fork);
		unless($pid)
		{
			open(IN,$bioms_thread[$c]);
                        my $out = $bioms_thread[$c] .'.l75';
                        open(OUT,">$out");
                        my @block;
                        my $c = 0;
                        while(<IN>)
                        {
                               $c ++;
                               chomp(my $line = $_);
                               if($line=~/^\@A0/)
                               {
                                       @block = ();
                                       push(@block,$line);
                               }
                               else
                               {
                                       push(@block,$line);
                                       if($c == 4)
                                       {
                                               if(length($block[1]) >= 75)
                                               {
                                                       my $nline = join("\n",@block);
                                                       print OUT "$nline\n";
                                               }
                                               $c = 0;
                                       }
                               }
                        }
                        close IN;
                        close OUT;
                        exit(0);
                }
	}
	foreach (@pid) 
        {
	        waitpid($_,0);
        }
}
#################################
sub change_headers
{
	(my @data) = @_;
	my $in = $data[0];
	my $out = $data[1];
	my $pre = $data[2];
	my $rem = $data[3];
	open(OUT,">$out");
	open(FAS,$in);
	while(<FAS>)
	{
		if($_ =~ /^>(.*)/)
		{
			chomp(my $name = $1);
			$name =~ s/$rem//;
			my $nname = $pre.'_'.$name;
			print OUT ">$nname\n";
		}
		else
		{
			print OUT "$_";
		}

	}
	close FAS;
	close OUT;	
}   
#################################
#################################
sub clean_that_fas
{
        (my @data) = @_;
        my $in = $data[0];
        my $out = $data[1];
        open(OUT,">$out");
        open(FAS,$in);
	while(<FAS>)
        {
                if($_ !~ /^>/)
                {
			chomp(my $seq = $_);
                        $seq =~ s/\*//g;
                        print OUT "$seq\n";
                }
                else
                {
			print OUT "$_";
		}

        }
        close FAS;
        close OUT;
}
#################################
#################################
#################################
#################################
sub prescreen
{
        (my @data) = @_;
        (my @bioms) = <$data[0]/*.split>;
        my @pid;
        my $db = '/media/yves/Xtra/MG-RAST_like/DEF/usefull_DBs/md5rRNA.C90.fasta';
        system "mkdir temp";
        my $num_threads = $data[1];
        my $nbtot = @bioms;
        my $calc = 0;
        for(my $j=0;$j<@bioms;$j += $num_threads)
        {
	        my $remain = $nbtot - $calc;
	        my $p = printf("%.2f", ($calc / $nbtot) * 100);	
	        print "$p % done\n";
	        if($remain > $num_threads)
	        {
		        my @bioms_thread;
		        for(my $k=$j;$k < $num_threads + $j;$k++)
		        {
			        push(@bioms_thread,$bioms[$k]);
		        }
		        my $c = -1;
		        for(my $i=0; $i < $num_threads; $i++)
		        {
			        $c ++;
			        $calc ++;
			        push(@pid,my $pid=fork);
			        unless($pid)
			        {
					(my $o) = ($bioms_thread[$c] =~ /\/([^\/]+)\.split/);
				        my $ficout = 'temp/'.$o.'.out';
				        system "./bin/blat $db $bioms_thread[$c] -out=blast8 -minIdentity=90  $ficout";
				        exit(0);
			        }
		        }
		        foreach (@pid) 
		        {
			        waitpid($_,0);
		        }
	        }
	        else
	        {
		        my @bioms_thread;
		        for(my $k=$j;$k < $nbtot;$k++)
		        {

			        push(@bioms_thread,$bioms[$k]);
		        }
		        my $c = -1;
		        for(my $i=0; $i < $remain; $i++)
		        {
			        $c ++;
			        push(@pid,my $pid=fork);
			        unless($pid)
			        {
				        (my $o) = ($bioms_thread[$c] =~ /\/([^\/]+)\.split/);
				        my $ficout = 'temp/'.$o.'.out';
				        system "./bin/blat $db $bioms_thread[$c] -out=blast8 -minIdentity=90 $ficout";
				        exit(0);
			        }
			
		        }
		        foreach (@pid) 
		        {
			        waitpid($_,0);
		        }
	        }
        }
}
###############################
###############################
###############################
###############################
sub blat_it_faa
{
        (my @data) = @_;
        system "mkdir temp";
        (my @bioms) = <md5nr_split2/*.faa>;
        my @pid;
        my $db = $data[0];
        my $num_threads = $data[1];
        my $nbtot = @bioms;
        my $calc = 0;
        for(my $j=0;$j<@bioms;$j += $num_threads)
        {
                my $remain = $nbtot - $calc;
                my $p = printf("%.2f", ($calc / $nbtot) * 100);	
                print "$p % done\n";
                if($remain > $num_threads)
                {
	                my @bioms_thread;
	                for(my $k=$j;$k < $num_threads + $j;$k++)
	                {
		                if(-s $bioms[$k])
		                {
		                        push(@bioms_thread,$bioms[$k]);
	                        }
	                }
	                my $c = -1;
	                for(my $i=0; $i < $num_threads; $i++)

	                {
		                $c ++;
		                $calc ++;
		                push(@pid,my $pid=fork);
		                unless($pid)
		                {
			                (my $o) = ($bioms_thread[$c] =~ /(\d+)\.faa/);
			                my $ficout = 'temp/'.$o.'.out';
			                #print "\n\n## $element vs $bioms_thread[$c] ##\n\n";
			                system "./bin/blat $bioms_thread[$c] $db -out=blast8 -minIdentity=60 -prot $ficout > blat_log";
			                exit(0);
		                }
	                }
	                foreach (@pid) 
	                {
		                waitpid($_,0);
	                }
                }
                else
                {
	                my @bioms_thread;
	                for(my $k=$j;$k < $nbtot;$k++)
	                {
		                if(-s $bioms[$k])
		                {
		                        push(@bioms_thread,$bioms[$k]);
	                        }
	                }
	                my $c = -1;
	                for(my $i=0; $i < $remain; $i++)
	                {
		                $c ++;
		                push(@pid,my $pid=fork);
		                unless($pid)
		                {
			                (my $o) = ($bioms_thread[$c] =~ /(\d+)\.faa/);
			                my $ficout = 'temp/'.$o.'.out';
			                #print "\n\n## $element vs $bioms_thread[$c] ##\n\n";
			                system "./bin/blat $bioms_thread[$c] $db -out=blast8 -minIdentity=60 -prot $ficout > blat_log";
			                exit(0);
		                }
		
	                }
	                foreach (@pid) 
	                {
		                waitpid($_,0);
	                }
                }
        }
}
#################################
#################################
sub blat_it_faa_split
{
        (my @data) = @_;
        system "mkdir temp";
        (my @bioms) = <md5nr_split/*.faa>;
        my @pid;
        (my @db) = <$data[0]/*.split>;
        my $num_threads = $data[1];
        my $nbtot = @bioms * @db;
        my $calc = 0;
        foreach my $element(@db)
        {
                for(my $j=0;$j<@bioms;$j += $num_threads)
                {
	                my $remain = $nbtot - $calc;
	                my $p = printf("%.2f", ($calc / $nbtot) * 100);	
	                print "$p % done\n";
	                if($remain > $num_threads)
	                {
		                my @bioms_thread;
		                for(my $k=$j;$k < $num_threads + $j;$k++)
		                {
			                if(-s $bioms[$k])
			                {
			                        push(@bioms_thread,$bioms[$k]);
		                        }
		                }
		                my $c = -1;
		                for(my $i=0; $i < $num_threads; $i++)
		                {
			                $c ++;
			                $calc ++;
			                push(@pid,my $pid=fork);
			                unless($pid)
			                {
				                (my $o) = ($bioms_thread[$c] =~ /(\d+)\.faa/);
				                my $ficout = 'temp/'.$o.'.out';
				                #print "\n\n## $element vs $bioms_thread[$c] ##\n\n";
				                system "./bin/blat $bioms_thread[$c] $element -out=blast8 -minIdentity=60 -prot $ficout > blat_log";
				                exit(0);
			                }
		                }
		                foreach (@pid) 
		                {
			                waitpid($_,0);
		                }
	                }
	                else
	                {
		                my @bioms_thread;
		                for(my $k=$j;$k < $nbtot;$k++)
		                {
			                if(-s $bioms[$k])
			                {
			                        push(@bioms_thread,$bioms[$k]);
		                        }
		                }
		                my $c = -1;
		                for(my $i=0; $i < $remain; $i++)
		                {
			                $c ++;
			                push(@pid,my $pid=fork);
			                unless($pid)
			                {
				                (my $o) = ($bioms_thread[$c] =~ /(\d+)\.faa/);
				                my $ficout = 'temp/'.$o.'.out';
				                #print "\n\n## $element vs $bioms_thread[$c] ##\n\n";
				                system "./bin/blat $bioms_thread[$c] $element -out=blast8 -minIdentity=60 -prot $ficout > blat_log";
				                exit(0);
			                }
			
		                }
		                foreach (@pid) 
		                {
			                waitpid($_,0);
		                }
	                }
                }
        }
}
#################################
#################################
#################################
#################################
sub blat_it_fna_split2
{
        (my @data) = @_;
        system "mkdir temp";
        (my @bioms) = <md5rna_split/*.fna>;
        my @pid;
        (my @db) = <$data[0]/*.split>;
        my $num_threads = $data[1];
        my $nbtot = @bioms * @db;
        my $calc = 0;
        foreach my $element(@db)
        {
                for(my $j=0;$j<@bioms;$j += $num_threads)
                {
	                my $remain = $nbtot - $calc;
	                my $p = printf("%.2f", ($calc / $nbtot) * 100);	
	                print "$p % done\n";
	                if($remain > $num_threads)
	                {
		                my @bioms_thread;
		                for(my $k=$j;$k < $num_threads + $j;$k++)
		                {
			                if(-s $bioms[$k])
			                {
			                        push(@bioms_thread,$bioms[$k]);
		                        }
		                }
		                my $c = -1;
		                for(my $i=0; $i < $num_threads; $i++)
		                {
			                $c ++;
			                $calc ++;
			                push(@pid,my $pid=fork);
			                unless($pid)
			                {
						if($bioms_thread[$c])
						{
							#print "\t\t\t $element VS $bioms_thread[$c] \n";
						        (my $o) = ($bioms_thread[$c] =~ /(\d+)\.fna/);
						        my $ficout = 'temp/'.$calc.'.out';
						        system "./bin/blat $bioms_thread[$c] $element -out=blast8 -minIdentity=97 $ficout";
						        exit(0);
					        }
			                }
		                }
		                foreach (@pid) 
		                {
			                waitpid($_,0);
		                }
	                }
	                else
	                {
		                my @bioms_thread;
		                for(my $k=$j;$k < $nbtot;$k++)

		                {
			                if(-s $bioms[$k])
			                {
			                        push(@bioms_thread,$bioms[$k]);
		                        }
		                }
		                my $c = -1;
		                for(my $i=0; $i < $remain; $i++)
		                {
			                $c ++;
			                push(@pid,my $pid=fork);
			                unless($pid)
			                {
						if($bioms_thread[$c])
						{
							#print "\t\t\t $element VS $bioms_thread[$c] \n";
						        (my $o) = ($bioms_thread[$c] =~ /(\d+)\.fna/);
						        my $ficout = 'temp/'.$o.'.out';
						        #print "\n\n## $element vs $bioms_thread[$c] ##\n\n";
						        system "./bin/blat $bioms_thread[$c] $element -out=blast8 -minIdentity=97 $ficout > rRNA_log";
						        exit(0);
					        }
			                }
			
		                }
		                foreach (@pid) 
		                {
			                waitpid($_,0);
		                }
	                }
                }
        }
}
#################################
#################################
sub blat_it_fna
{
        (my @data) = @_;
        system "mkdir temp2";
        (my @bioms) = <md5rna_split/*.fna>;
        my @pid;
        my $db = $data[0];
        my $num_threads = $data[1];
        my $nbtot = @bioms;
        my $calc = 0;
        for(my $j=0;$j<@bioms;$j += $num_threads)
        {
	        my $remain = $nbtot - $calc;
	        my $p = printf("%.2f", ($calc / $nbtot) * 100);	
	        print "$p % done\n";
	        if($remain > $num_threads)
	        {
		        my @bioms_thread;
		        for(my $k=$j;$k < $num_threads + $j;$k++)
		        {
			        push(@bioms_thread,$bioms[$k]);
		        }
		        my $c = -1;
		        for(my $i=0; $i < $num_threads; $i++)
		        {
			        $c ++;
			        $calc ++;
			        push(@pid,my $pid=fork);
			        unless($pid)
			        {
				        (my $o) = ($bioms_thread[$c] =~ /\/([^\/]+)\.fna/);
				        my $ficout = 'temp2/'.$o.'.out';
				        system "./bin/blat $db $bioms_thread[$c] -out=blast8 -minIdentity=97 $ficout";
				        exit(0);
			        }
		        }
		        foreach (@pid) 
		        {
			        waitpid($_,0);
		        }
	        }
	        else
	        {
		        my @bioms_thread;
		        for(my $k=$j;$k < $nbtot;$k++)
		        {

			        push(@bioms_thread,$bioms[$k]);
		        }
		        my $c = -1;
		        for(my $i=0; $i < $remain; $i++)
		        {
			        $c ++;
			        push(@pid,my $pid=fork);
			        unless($pid)
			        {
				        (my $o) = ($bioms_thread[$c] =~ /\/([^\/]+)\.fna/);
				        my $ficout = 'temp2/'.$o.'.out';
				        system "./bin/blat $db $bioms_thread[$c] -out=blast8 -minIdentity=97 $ficout";
				        exit(0);
			        }
			
		        }
		        foreach (@pid) 
		        {
			        waitpid($_,0);
		        }
	        }
        }
}
#################################
#################################
sub subsample_fasta
{
        (my @data) = @_;
        #(my @subsamples) = (0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1);
        #(my @subsamples) = (0.1,0.2,0.3,0.4,0.5,0.75);
        #(my @subsamples) = (0.001,0.01,0.05,1);
	(my @subsamples) = (5000000);
        my $ref = $data[0];
        my $pre = $data[1];
        my $in = $data[2];
	my $name;
	my %fasta;
	my $tot = 0;
	open(IN,$in);
	while(<IN>)
        {
                if($_ =~ /^>/)
                {
                        chomp($name = $_);
                        $tot ++;
                }
                else
                {
                        chomp($fasta{$name} .= $_);
                }
        }
        close IN;
        foreach my $percent(@subsamples)
        {
                for(my $repet = 1; $repet <= 1; $repet ++)
                {

                     print "##### Subsampling $percent Repetition $repet ####\n\n";
                     my $out = $pre.'/Assembly_IDBA/'.$pre.'_'.$percent.'_'.$repet.'.fasta';
                     open(OUT,">$out");
                     #my $frac = int($percent * $tot);
		     my $frac = $percent;
                     print "\n\n\t\t\tshould retrieve $frac seqs\n\n";
                     my $count = 0;
                     foreach my $k(keys(%fasta))
                     {
                             if($count <= $frac)
                             {
                                my $ac = rand(1);
                                if($ac <= 0.5)
                                {
                                	print OUT "$k\n$fasta{$k}\n";
                                	$count ++
                               	}
                             }
                             else
                             {
                                last;
                             }
                     }
                     close OUT;
                }
        }
        %fasta = ();
        my $out = $pre.'/Assembly_IDBA/'.$pre.'_1.fasta';
        system "cp $in $out";
}
#################################
#################################
sub subsample_fasta2
{
        (my @data) = @_;
        my $ref = $data[0];
        my $pre = $data[1];
        my $in = $data[2];
        my $percent = $data[3];
        my $repet = $data[4];
	my $name;
	my %fasta;
	my $tot = 0;
	open(IN,$in);
	while(<IN>)
        {
                if($_ =~ /^>/)
                {
                        chomp($name = $_);
                        $tot ++;
                }
                else
                {
                        chomp($fasta{$name} .= $_);
                }
        }
        close IN;
        my $out = $pre.'/Assembly_IDBA/'.$ref.'_'.$percent.'_'.$repet.'.fasta';
        open(OUT,">$out");
        my $frac = $percent * $tot;
        my $count = 0;
        foreach my $k(keys(%fasta))
        {
                $count ++;
                if($count <= $frac)
                {
                       print OUT "$k\n$fasta{$k}\n";
                }
                else
                {
                       last;
                }
       }
       close OUT;
}
#################################
#################################
#################################
#################################
sub subsample_fastq
{
        (my @data) = @_;
        my $in1 = $data[0];
        my $in2 = $data[1];
        my $per = $data[2];
	my @pid;
        my $nbtot = 2;
	my @bioms_thread;
	push(@bioms_thread,$in1,$in2);
	my $c = -1;
	for(my $i=0; $i < 2; $i++)
	{
		$c ++;
		push(@pid,my $pid=fork);
		unless($pid)
		{
			my $out = $bioms_thread[$c].'_ss_'.$per;
			open(IN,$bioms_thread[$c]);
	                open(OUT, ">$out");
	                my $cl = 0;
	                my $count = 0;
	                while(<IN>)
                        {
                                $count ++;
                        }
                        close IN;
                        open(IN,$bioms_thread[$c]);
                        my $subl = $count * $per;
                        my %fq;
                        $cl = 0;
                        while(<IN>)
                        {
                                $cl ++;
                                my $line = $_;
                                if($cl < $subl)
                                {
                                       print OUT "$line";
                                }
                                else
                                {
                                       last;
                                }
                        }
                        close IN;
                        close OUT;	
			exit(0);
		}
	}
	foreach (@pid) 
        {
	        waitpid($_,0);
        }
}
#################################
#################################
sub split_da_file_C90
{
        (my @data) = @_;
        my $in = $data[0];
        my $nb = $data[1];
        print "\n\n#### Splitting $in / $nb seqs per file ####\n\n";
        system "mkdir $data[2]";
	my $t = 0;
	my $cfile = 1;
	my $fileout = $data[2].'/'.$cfile.'.split';
	open(IN,$in);
	open(OUT, ">$fileout");
	while(<IN>)
	{		
		if($_ =~ /^>/)
		{
		        $t ++;
		        if($t >= $nb)
		        {
		                $t = 0;
		                close OUT;
			        $cfile ++;
			        $fileout =  $data[2].'/'.$cfile.'.split';
			        open(OUT, ">$fileout");
		        }        
		}	        
		print OUT "$_";
	}
	close IN;
	close OUT;
}
#################################
#################################
sub split_da_file
{
        (my @data) = @_;
        my $in = $data[0];
        my $nb = $data[1];
        system "mkdir $data[2]";
        print "\n\n#### Splitting $in / $nb seqs per file at $data[2] path####\n\n";
	my $t = 0;
	my $cfile = 1;
	my $fileout = $data[2].'/'.$cfile.'.split';
	open(IN,$in);
	open(OUT, ">$fileout");
	while(<IN>)
	{		
		
		if($_ =~ /^>/)
		{
		        $t ++;
		        if( ($t % $nb) == 0 )
		        {
			        close OUT;
			        $cfile ++;
			        $fileout = $data[2].'/'.$cfile.'.split';
			        open(OUT, ">$fileout");
		        }	        
		        print OUT "$_";
	        }
	        else
	        {
	                print OUT "$_";
                }
		
	}
	close IN;
	close OUT;
}
#################################
#################################
sub split_da_file_and_convert
{
        (my @data) = @_;
        my $in = $data[0];
        my $nb = $data[1];
        system "mkdir $data[2]";
        print "\n\n#### Splitting $in / $nb seqs per file at $data[2] path####\n\n";
	my $t = 0;
	my $cfile = 1;
	my $fileout = $data[2].'/'.$cfile.'.split';
	open(IN,$in);
	open(OUT, ">$fileout");
	while(<IN>)
	{		
		
		if($_ =~ /^>/)
		{
		        $t ++;
		        chomp(my $line = $_);
		        if( ($t % $nb) == 0 )
		        {
			        close OUT;
			        $cfile ++;
			        $fileout = $data[2].'/'.$cfile.'.split';
			        open(OUT, ">$fileout");
		        }
		        $line =~ s/\s+/#/;        
		        print OUT "$line\n";
	        }
	        else
	        {
	                print OUT "$_";
                }
		
	}
	close IN;
	close OUT;
}
#################################
#################################
sub abundance
{
        (my @data) = @_;
	open(IN,$data[0]);
	my %temp;
	my %refs;
	my $name;
	while(<IN>)
	{
		if($_ =~ /^>(.*)/)
		{
			chomp($name = $1);
		}
		else
		{
			$temp{$name} ++;
			if($_ =~ />(.*)\.\.\. \*/)
			{
				$refs{$name} = $1;
			}
		}	
	}
	close IN;
	(my @s) = sort{$a cmp $b} keys(%temp);
	my %abd;
	foreach my $e(@s)
	{
		$abd{$refs{$e}} = $temp{$e};
	}
	return(%abd);
	@s=();
        %abd=();
        %temp=();
}
#################################
#################################
sub make_directories
{
        (my @data) = @_;
        my $pre = 'BeTAE_'.$data[0];
        if(not -d $pre)
        {
                system "mkdir $pre";
        }
        my $SolexaQC = $pre.'/SolexaQC';
        if(not -d $SolexaQC)
        {
                system "mkdir $SolexaQC";
        }
        my $IDBA = $pre.'/Assembly_IDBA';
        if(not -d $IDBA)
        {
                system "mkdir $IDBA";
        }
        my $PROTS = $pre.'/Gene_Prediction_FragGeneScan';
        if(not -d $PROTS)
        {
                system "mkdir $PROTS";
        }
        my $SIM = $pre.'/Similarity_search';
        if(not -d $SIM)
        {
                system "mkdir $SIM";
        }
        my $ANNOTATION = $pre.'/Functions';
        if(not -d $ANNOTATION)
        {
                system "mkdir $ANNOTATION";
                my $pat = $ANNOTATION.'/COGs';
                system "mkdir $pat";
                $pat = $ANNOTATION.'/KOs';
                system "mkdir $pat";
                $pat = $ANNOTATION.'/SUBSYSTEMS';
                system "mkdir $pat";
        }
        my $TAXONOMY = $pre.'/Taxonomy';
        if(not -d $TAXONOMY)
        {
                system "mkdir $TAXONOMY";
                my $pat = $TAXONOMY.'/LCA';
                system "mkdir $pat";
                $pat = $TAXONOMY.'/Absolute';
                system "mkdir $pat";
        }
        my $rRNA = $pre.'/rRNA';
        if(not -d $rRNA)
        {
                my $pat = $rRNA;
                system "mkdir $pat";
        }
}
#################################
#################################
sub stats_fasta
{
        (my @data) = @_;
        my $in = $data[0];
        open(IN,$in);
        my $c = 0;
        my $l = 0;
        while(<IN>)
        {
                if($_ =~ /^>/)
                {
                        $c ++;
                }
                else
                {
                        chomp(my $line = $_);
                        $line =~ s/\s//g;
                        $l += length($line);
                }
        }
        close IN;
        return($c,$l);
}
#################################
#################################
sub stats_sim
{
        (my @data) = @_;
        my $in = $data[0];
        open(IN,$in);
        my %res;
        while(<IN>)
        {
                (my @tab) = split("\t",$_);
                $res{$tab[0]} ++;
        }
        close IN;
        (my $c) = keys(%res);
        return($c);
}
#################################
#################################
sub stats_fastq
{
        (my @data) = @_;
        open(IN,$data[0]);
        my @block;
        my $cline = 0;
        my $c = 0;
        my $l = 0;
        while(<IN>)
        {
                $cline ++;
                chomp(my $line = $_);
                if($line=~/^\@A0/)
                {
                        @block = ();
                        push(@block,$line);
                        $c ++;
                }
                else
                {
                        push(@block,$line);
                        if($cline == 4)
                        {
                                $block[1] =~ s/\s//g;
                                $l += length($block[1]);
                                $cline = 0;
                        }
                 }
        }
        close IN;
        return($c,$l);
}
#################################
#################################     




#### tips ####

	### KRAKEN ###
	# /media/yves/Xtra/MG-RAST_like/DEF/usefull_DBs/ftp.metagenomics.anl.gov/FINAL$ ../../../kraken2-2.0.7-beta/kraken2-build --download-taxonomy --threads 28 --db rRNA_profiles.db
	# kraken2-2.0.7-beta/kraken2-build --add-to-library md5rRNA.kraken --db rRNA_profiles.db --threads 28
	# kraken2-build --build --db rRNA_profiles.db
