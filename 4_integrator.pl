use strict;
use warnings;
use Getopt::Long;

my $gene = "";
my $GeMoMa_file = "";
my $PASA_file = "";
my $Braker_file = "";
my $genome = "";
my $mRNA = "";
my $output = "";
my $help = "";
my $overlap_ratio = 0.7;
my $number = 1;
my %genome = ();

my %gene_PASA = ();
my %gene_GeMoMa = ();
my %gene_Braker = ();

my %gene = ();
my %coding = ();
my %strand = ();
my $usage = "Usage: perl 4_integrator.pl -g <genome_info> --est <PASA_file> --homology <GeMoMa_file> --denovo <Braker_file> -o <output>\n";
GetOptions ("g=s" => \$genome,
			"est=s" => \$PASA_file,
			"homology=s" => \$GeMoMa_file,
			"denovo=s" => \$Braker_file,
			"r=f" => \$overlap_ratio,
			"o=s" => \$output,
			"h" => \$help) or die(get_help($usage));
if ($help) {get_help(); die()}
if (!$GeMoMa_file or !$genome or !$PASA_file or !$Braker_file or !$output) {get_help(); die($usage)}

open (INFO,"<$genome") or die("Please enter the genomic file: -g <genome_info>\n");
open (PASA,"<$PASA_file") or die("Please enter the PASA input file: --homology <PASA_gff3>\n");
open (GEMOMA,"<$GeMoMa_file") or die("Please enter the GeMoMa input file: --est <GeMoMa_gff3>\n");
open (BRAKER,"<$Braker_file") or die("Please enter the Braker input file: --denovo <Braker_gff3>\n");;
open (OUT,">$output") or die("Please specify the output file: -o <output>\n");

print "Loading genome info...";
while (<INFO>) {
	chomp;
	my @tmp = split /\t/,$_;
	$genome{$tmp[0]} = $tmp[1];
}
print " there are ",scalar keys %genome," contigs.\n";

print "Loading PASA-annotated mRNAs... ";
my $PASA_number = 0;
while (<PASA>) {
	chomp;
	if (/^#/ or !/\S/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
	if ($tmp[2] eq "gene" and $tmp[8] =~ /ID=(\S+)/) {
		$gene = $1;
		$number = 1;
		$mRNA = "";
	}
	if ($tmp[2] eq "mRNA" and $tmp[8] =~ /ID=(\S+?);/) {
		$PASA_number++;
		$mRNA = $1;
		push @{$gene_PASA{$chr}{$gene}}, $mRNA;
		$strand{$mRNA} = $tmp[6];
	}
	elsif ($tmp[2] eq "CDS" and $tmp[8] =~ /Parent=(\S+)/) {
		push @{$coding{$mRNA}}, $_;
	}
}
print " there are $PASA_number PASA mRNAs.\n";

print "Loading GeMoMa-annotated mRNAs... ";
my $GeMoMa_number = 0;
while (<GEMOMA>) {
	chomp;
	if (/^#/ or !/\S/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
	if ($tmp[2] eq "gene" and $tmp[8] =~ /ID=(\S+)/) {
		$gene = $1;
		$number = 1;
		$mRNA = "";
	}
	if ($tmp[2] eq "mRNA" and $tmp[8] =~ /ID=(\S+?);/) {
		$GeMoMa_number++;
		$mRNA = $1;
		push @{$gene_GeMoMa{$chr}{$gene}}, $mRNA;
		$strand{$mRNA} = $tmp[6];
	}
	elsif ($tmp[2] eq "CDS" and $tmp[8] =~ /Parent=(\S+)/) {
		push @{$coding{$mRNA}}, $_;
	}
}
print " there are $GeMoMa_number GeMoMa mRNAs.\n";

print "Loading Braker-annotated mRNAs... ";
my $Braker_number = 0;
while (<BRAKER>) {
	chomp;
	if (/^#/ or !/\S/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
	if ($tmp[2] eq "gene" and $tmp[8] =~ /ID=(\S+)/) {
		$gene = $1;
		$number = 1;
		$mRNA = "";
	}
	if ($tmp[2] eq "mRNA" and $tmp[8] =~ /ID=(\S+?);/) {
		$Braker_number++;
		$mRNA = $1;
		push @{$gene_Braker{$chr}{$gene}}, $mRNA;
		$strand{$mRNA} = $tmp[6];
	}
	elsif ($tmp[2] eq "CDS" and $tmp[8] =~ /Parent=(\S+)/) {
		push @{$coding{$mRNA}}, $_;
	}
}
print " there are $Braker_number Braker mRNAs.\n";
$number = 1;

foreach my $chr (sort keys %genome) {
	print "Processing $chr... ";
	my $anno = "0" x $genome{$chr};
	my %anno = ();
	my %pass = ();
	my %region = ();
	my %length = ();
	my %type = ();
	# map the PASA annotated mRNAs onto the genomic sequence
	foreach my $gene (sort keys %{$gene_PASA{$chr}}) {
		foreach my $mRNA (@{$gene_PASA{$chr}{$gene}}) {
			my ($begin,$end,$length) = get_region($mRNA,\@{$coding{$mRNA}});
			$region{$mRNA} = $begin."\t".$end;
			$length{$mRNA} = $length;
			foreach my $CDS (@{$coding{$mRNA}}) {
				my @tmp = split /\t/,$CDS;
				substr($anno,$tmp[3],$tmp[4]-$tmp[3]+1) = "2" x ($tmp[4]-$tmp[3]+1);
				for (my $i = $tmp[3]; $i < $tmp[4]; $i++) {$anno{$i}{$mRNA} = 1}
			}
			$pass{$mRNA} = $begin;
			$type{$mRNA} = "PASA";
		}
	}
	# map the GeMoMa annotated mRNAs onto the genomic sequence
	foreach my $gene (sort keys %{$gene_GeMoMa{$chr}}) {
		foreach my $mRNA (@{$gene_GeMoMa{$chr}{$gene}}) {
			my ($begin,$end,$length) = get_region($mRNA,\@{$coding{$mRNA}});
			$region{$mRNA} = $begin."\t".$end;
			$length{$mRNA} = $length;

			# calculate the overlapping ratio between the current mRNA and PASA annotated ones
			my $string = "";
			foreach my $CDS (@{$coding{$mRNA}}) {
				my @tmp = split /\t/,$CDS;
				$string .= substr($anno,$tmp[3],$tmp[4]-$tmp[3]+1);
			}
			my @cover = ();
			@cover = $string =~ /2/g;

			# if the overlapping ratio < 10%, then regard the current mRNA as confident, and map this mRNA onto the genomic seqeunce
			if (scalar(@cover)/length($string) < 0.1) {
				foreach my $CDS (@{$coding{$mRNA}}) {
					my @tmp = split /\t/,$CDS;
					substr($anno,$tmp[3],$tmp[4]-$tmp[3]+1) = "2" x ($tmp[4]-$tmp[3]+1);
					for (my $i = $tmp[3]; $i < $tmp[4]; $i++) {$anno{$i}{$mRNA} = 1}
				}
				$pass{$mRNA} = $begin;
				$type{$mRNA} = "GeMoMa";
			}
			# if the overlapping ratio >= 10%, then compare current mRNA with newly annotated ones
			else {
				my %old = ();
				my %old_exon = ();
				my $alter = 0;
				# get the list of overlapped mRNA
				foreach my $CDS (@{$coding{$mRNA}}) {
					my @tmp = split /\t/,$CDS;
					my %count = ();
					for (my $i = $tmp[3]; $i < $tmp[4]; $i++) {
						if (exists $anno{$i}) {
							foreach my $old (sort keys %{$anno{$i}}) {
								if ($pass{$old} ne "fail") {
									$count{$old}++;
									$old{$old}++;
								}
							}
						}
					}
					# for each exon of the current mRNA, if previously annotated mRNAs cover >50% sequence, consider a exon-exon synteny between them 
					foreach my $old (sort keys %count) {
						if ($count{$old} / ($tmp[4]-$tmp[3]+1) > 0.5) {$old_exon{$old}++} else {$old_exon{$old} += 0}
					}					
				}
				# for each overlapped previously annotated mRNAs, if >70% sequence or >50%(>2) exons are covered by current mRNA, use the current mRNA to 
				# replace it and retreive the previously mapped genomic region
				foreach my $old (sort keys %old) {
					if ($length{$old} <= $length{$mRNA}) {
						if ($old{$old}/$length{$old} >= $overlap_ratio or ($old_exon{$old}/scalar(@{$coding{$old}}) >= 0.5 and $old_exon{$old} >= 2)) {
							foreach my $CDS (@{$coding{$old}}) {
								my @tmp = split /\t/,$CDS;
								substr($anno,$tmp[3],$tmp[4]-$tmp[3]+1) = "0" x ($tmp[4]-$tmp[3]+1);
							}
							$pass{$old} = "fail";
							$alter++;
						}
					}
				}

				# re-map the current mRNA onto the genomic seqeunce after removing, if any, the "failed" mRNAs 
				$string = "";
				foreach my $CDS (@{$coding{$mRNA}}) {
					my @tmp = split /\t/,$CDS;
					$string .= substr($anno,$tmp[3],$tmp[4]-$tmp[3]+1);
				}
				my @cover = ();
				@cover = $string =~ /2/g;

				# if >70% sequence of the current mRNA are cover by previously annotated ones, discard this mRNA
				if (scalar(@cover)/length($string) >= $overlap_ratio) {next}

				# map the current mRNA onto the genomic sequence if it conforms to the criterion above
				foreach my $CDS (@{$coding{$mRNA}}) {
					my @tmp = split /\t/,$CDS;
					substr($anno,$tmp[3],$tmp[4]-$tmp[3]+1) = "2" x ($tmp[4]-$tmp[3]+1);
					for (my $i = $tmp[3]; $i < $tmp[4]; $i++) {$anno{$i}{$mRNA} = 1}
				}
				$pass{$mRNA} = $begin;
				$type{$mRNA} = "GeMoMa";
			}
		}
	}
	# release the memory
	%anno = ();

	# map the Braker annotated mRNAs onto the genomic sequence
	foreach my $gene (sort keys %{$gene_Braker{$chr}}) {
		foreach my $mRNA (@{$gene_Braker{$chr}{$gene}}) {
			my ($begin,$end,$length) = get_region($mRNA,\@{$coding{$mRNA}});
			$region{$mRNA} = $begin."\t".$end;
			$length{$mRNA} = $length;

			# calculate the overlapping ratio between the current mRNA and PASA annotated ones
			my $string = "";
			foreach my $CDS (@{$coding{$mRNA}}) {
				my @tmp = split /\t/,$CDS;
				$string .= substr($anno,$tmp[3],$tmp[4]-$tmp[3]+1);
			}
			my @cover = ();
			@cover = $string =~ /2/g;

			# if the overlapping region < 100 bp, then keep it
			if (scalar @cover < 100) {
				$pass{$mRNA} = $begin;
				$type{$mRNA} = "Braker";
			}
		}
	}

	#prepare output
	my @final = grep {$pass{$_} ne "fail"} keys %pass;
	print "done!\n";
	print "Writing gff3 file for $chr... ";
	my $output_number = 0;
	foreach my $mRNA (sort {$pass{$a} <=> $pass{$b}} @final) {
		$output_number++;
		my $type = $type{$mRNA};
		print OUT "$chr\t$type\tgene\t",$region{$mRNA},"\t.\t",$strand{$mRNA},"\t.\tID=gene_$number\n";
		print OUT "$chr\t$type\tmRNA\t",$region{$mRNA},"\t.\t",$strand{$mRNA},"\t.\tID=$mRNA;Parent=gene_$number\n";
		print OUT join("\n",@{$coding{$mRNA}}),"\n\n";
		$number++;
	}
	print "a total of $output_number mRNAs are kept!\n";
}

sub get_region {
	my ($i,$j) = @_;
	my @tmp = @$j;
	my ($begin,$end,$length) = (0,0,0);
	foreach my $cds (@tmp) {
		my @tmp2 = split /\t/,$cds;
		$length += $tmp2[4]-$tmp2[3]+1;
	}
	if ($strand{$i} eq "+") {
		my @tmp2 = split /\t/,$tmp[0];
		$begin = $tmp2[3];
		@tmp2 = split /\t/,$tmp[-1];
		$end = $tmp2[4];
	}
	elsif ($strand{$i} eq "-") {
		my @tmp2 = split /\t/,$tmp[-1];
		$begin = $tmp2[3];
		@tmp2 = split /\t/,$tmp[0];
		$end = $tmp2[4];
	}
	return ($begin,$end,$length);
}

sub get_help {
print <<END_OF_TEXT;

    ###########################
    #        Runxian Yu       #
    #  yurunxian\@ibcas.ac.cn  #
    ###########################

    This script is designed for integrating PASA (transcript), GeMoMa (homology) and Braker (de novo) gene
    preditions onto the gemonic sequence.

    Please use the output gff3 of 1_PASA_parser.pl, 2_GeMoMa_parser.pl and 3_Braker_parser.pl!!!

    Usage:
    -g <genome_info>      a text file containing the name and length for each contig/scaffold/chromosome
                          example:    chr1    1000000
                                      chr2    1000000
                                      chr3    1000000
    --est      <PASA_gff3>
                          PASA gff3 output of 1_PASA_parser.pl
    --homology <GeMoMa_gff3>
                          GeMoMa gff3 output of 2_GeMoMa_parser.pl
    --denovo   <Braker_gff3>
                          Braker gff3 output of 3_Braker_parser.pl
    -r <0-1>              the overlapping ratio used to compare mRNAs. Default: 0.7
    -o <output>           output file
    -h                    print this page

END_OF_TEXT
}
