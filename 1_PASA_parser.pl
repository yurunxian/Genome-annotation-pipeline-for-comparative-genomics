use strict;
use warnings;
use Getopt::Long;

my $gene = "";
my $PASA_file = "";
my $genome = "";
my $mRNA = "";
my $output = "";
my $help = "";
my $overlap_ratio = 0.7;
my $remove_single_exon = "";
my $number = 1;
my %genome = ();
my %TE = ();
my %gene = ();
my %coding = ();
my %strand = ();
my $usage = "Usage: perl 1_PASA_parser.pl -g <genome_info> -i <PASA_file> -o <output>\n";
GetOptions ("i=s" => \$PASA_file,
			"g=s" => \$genome,
			"r=f" => \$overlap_ratio,
			"o=s" => \$output,
			"h" => \$help,
			"remove_single_exon" => \$remove_single_exon) or die(get_help());
if ($help) {get_help(); die()}
if ($remove_single_exon) {$remove_single_exon = 1}
if (!$PASA_file or !$genome or !$output) {get_help(); die($usage)}

open (INFO,"<$genome") or die($usage);
open (PASA,"<$PASA_file") or die($usage);
open (OUT,">$output") or die($usage);

print "Loading genome info...";
while (<INFO>) {
	chomp;
	my @tmp = split /\s+/,$_;
	$genome{$tmp[0]} = $tmp[1];
}
print " there are ",scalar keys %genome," contigs.\n";

print "Loading PASA file... ";
my $PASA_number = 0;
while (<PASA>) {
	chomp;
	if (/^#/ or !/\S/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
	if ($tmp[2] eq "gene" and $tmp[8] =~ /ID=(\S+?);/) {
		$gene = $1;
	}
	elsif ($tmp[2] eq "mRNA" and $tmp[8] =~ /ID=(\S+?);/) {
		$PASA_number++;
		$mRNA = $1;
		push @{$gene{$chr}{$gene}}, $mRNA;
		$strand{$mRNA} = $tmp[6];
	}
	elsif ($tmp[2] eq "CDS" and $tmp[8] =~ /Parent=(\S+)/) {
		$mRNA = $1;
		push @{$coding{$mRNA}}, $_;
	}
}
print " a total of $PASA_number mRNAs are imported!\n";

foreach my $chr (sort keys %genome) {
	print "Processing $chr... ";
	my $anno = "0" x $genome{$chr};
	my $last = "";
	my %anno = ();
	my %pass = ();
	my %region = ();
	my %length = ();
	foreach my $gene (sort keys %{$gene{$chr}}) {
		foreach my $mRNA (@{$gene{$chr}{$gene}}) {
			my ($begin,$end,$length) = get_region($mRNA,\@{$coding{$mRNA}});
			$region{$mRNA} = $begin."\t".$end;
			$length{$mRNA} = $length;

			# remove any genes longer than 150kb and single-exon genes shorter than 300bp
			if ($end-$begin+1 > 150000) {next}
			if (scalar @{$coding{$mRNA}} eq 1 and $length{$mRNA} <= 300) {next}

			# calculate the overlapping ratio between the current mRNA and previously annotated ones
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
				$last = $mRNA;
			}
			# if the overlapping ratio >= 10%, then compare current mRNA with previously annotated ones
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

				# re-map the current mRNA onto the genomic seqeunce after clean, if any, the "failed" mRNAs 
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
				$last = $mRNA;
			}
		}
	}
	# release the memory
	%anno = ();
	my @final = grep {$pass{$_} ne "fail"} keys %pass;
	print "done!\n";
	
	# remove single-exon genes nested within multiple-exon ones
	if ($remove_single_exon eq 1) {
		my $remove_number = 0;
		print "Removing single-exon genes for $chr...";
		foreach my $mRNA (sort {$pass{$a} <=> $pass{$b}} @final) {
			my ($begin,$end) = split /\t/,$region{$mRNA};
			for (my $i = $begin; $i < $end; $i++) {$anno{$i}{$mRNA} = 1}
		}
		for (my $i = 1; $i < $genome{$chr}; $i += 10000) {
			my %overlap = ();
			for (my $j = $i; $j < $i+20000; $j++) {
				if (exists $anno{$j} and scalar keys %{$anno{$j}} > 1) {
					my @symbol = ();
					foreach my $k (sort keys %{$anno{$j}}) {if ($pass{$k} ne "fail") {push @symbol, $k}}
					if (scalar @symbol >= 2) {$overlap{join(",",@symbol)}++}
				}
			}
			my @overlap = grep {$overlap{$_} > 300} keys %overlap;
			if (scalar @overlap < 1) {next}
			foreach my $m (sort @overlap) {
				foreach my $n (split /,/,$m) {if (scalar @{$coding{$n}} eq 1) {$pass{$n} = "fail";$remove_number++}}
			}
		}
		print " a total of $remove_number single-exon mRNAs are removed!\n";
		# release the memory
		%anno = ();
		@final = grep {$pass{$_} ne "fail"} keys %pass;
	}

#prepare output
	print "Writing gff3 file for $chr... ";
	my $output_number = 0;
	foreach my $mRNA (sort {$pass{$a} <=> $pass{$b}} @final) {
		$output_number++;
		print OUT "$chr\ttransdecoder\tgene\t",$region{$mRNA},"\t.\t",$strand{$mRNA},"\t.\tID=PASA$number\n";
		print OUT "$chr\ttransdecoder\tmRNA\t",$region{$mRNA},"\t.\t",$strand{$mRNA},"\t.\tID=$mRNA;Parent=PASA$number\n";
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

    This script is designed for integrating PASA gene preditions onto the gemonic sequence.
    Please use the standard PASA output file!!!

    Usage:
    -g <genome_info>      a text file containing the name and length for each contig/scaffold/chromosome
                          example:    chr1    1000000
                                      chr2    1000000
                                      chr3    1000000
    -i <PASA_file>        the standard PASA output, e.g., compreh_init_build.fasta.transdecoder.genome.gff3
    -r <0-1>              the threshold of overlapping ratio used to merge two mRNA loci. Default: 0.7
    -o <output>           output file
    --remove_single_exon  remove single-exon mRNAs which are nested within multiple-exon ones
    -h                    print this page

END_OF_TEXT
}
