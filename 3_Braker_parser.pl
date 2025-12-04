use strict;
use warnings;
use Getopt::Long;

my $gene = "";
my $Braker_file = "";
my $genome = "";
my $mRNA = "";
my $output = "";
my $help = "";
my $overlap_ratio = 0.7;
my $number = 1;
my %genome = ();
my %gene = ();
my %coding = ();
my %strand = ();
my $usage = "Usage: perl 3_Braker_parser.pl -g <genome_info> -i <Braker_file> -o <output>\n";
GetOptions ("i=s" => \$Braker_file,
			"g=s" => \$genome,
			"r=f" => \$overlap_ratio,
			"o=s" => \$output,
			"h" => \$help) or die(get_help());
if ($help) {get_help(); die()}
if (!$Braker_file or !$genome or !$output) {get_help(); die($usage)}

open (INFO,"<$genome") or die($usage);
open (BRAKER,"<$Braker_file") or die($usage);
open (OUT,">$output") or die($usage);

print "Loading genome info...";
while (<INFO>) {
	chomp;
	my @tmp = split /\s+/,$_;
	$genome{$tmp[0]} = $tmp[1];
}
print " there are ",scalar keys %genome," contigs.\n";

print "Loading Braker file... ";
my $Braker_number = 0;
while (<BRAKER>) {
	chomp;
	if (/^#/ or !/\S/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
	if ($tmp[2] eq "gene" and $tmp[8] =~ /(\S+)/) {
		$gene = $1;
		$mRNA = "";
	}
	if ($tmp[2] eq "transcript" and $tmp[8] =~ /(\S+)/) {
		$Braker_number++;
		$mRNA = $1;
		push @{$gene{$chr}{$gene}}, $mRNA;
		$strand{$mRNA} = $tmp[6];
	}
	elsif ($tmp[2] eq "CDS") {
		$tmp[8] = "Parent=$mRNA";
		$tmp[1] = "Braker";
		my $string = join("\t",@tmp);
		if ($strand{$mRNA} eq "+") {push @{$coding{$mRNA}}, $string}
		elsif ($strand{$mRNA} eq "-") {unshift @{$coding{$mRNA}}, $string}
	}
}
print " a total of $Braker_number Braker mRNAs are imported!\n";

foreach my $chr (sort keys %genome) {
	print "Processing $chr... ";
	my $anno = "0" x $genome{$chr};
	my $last = "";
	my %anno = ();
	my %braker = ();
	my %region = ();
	my %length = ();
	foreach my $gene (sort keys %{$gene{$chr}}) {
		foreach my $mRNA (@{$gene{$chr}{$gene}}) {
			my ($begin,$end,$length) = get_region($mRNA,\@{$coding{$mRNA}});
			$region{$mRNA} = $begin."\t".$end;
			$length{$mRNA} = $length;

			# remove any genes longer than 150kb and single-exon genes
			if ($end-$begin+1 > 150000 or scalar @{$coding{$mRNA}} eq 1) {next}

			# calculate the overlapping ratio between the current mRNA and newly annotated ones
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
				$braker{$mRNA} = $begin;
				$last = $mRNA;
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
								if ($braker{$old} ne "fail") {
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
							$braker{$old} = "fail";
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
				$braker{$mRNA} = $begin;
				$last = $mRNA;
			}
		}
	}
	# release the memory
	%anno = ();
	my @final = grep {$braker{$_} ne "fail"} keys %braker;
	print "done!\n";

	#prepare output
	print "Writing gff3 file for $chr... ";
	my $output_number = 0;
	foreach my $mRNA (sort {$braker{$a} <=> $braker{$b}} @final) {
		$output_number++;
		print OUT "$chr\tBraker\tgene\t",$region{$mRNA},"\t.\t",$strand{$mRNA},"\t.\tID=Braker$number\n";
		print OUT "$chr\tBraker\tmRNA\t",$region{$mRNA},"\t.\t",$strand{$mRNA},"\t.\tID=$mRNA;Parent=Braker$number\n";
		print OUT join("\n",@{$coding{$mRNA}}),"\n\n";
		$number++;
	}
	print "a total of $output_number mRNAs are kept!\n";
}
#chr1	AUGUSTUS	gene	106802	107131	.	+	.	g_27000
#chr1	AUGUSTUS	transcript	106802	107131	1	+	.	anno1.g27818.t1
#chr1	AUGUSTUS	start_codon	106802	106804	.	+	0	transcript_id "anno1.g27818.t1"; gene_id "g_27000";
#chr1	AUGUSTUS	CDS	106802	107131	1	+	0	transcript_id "anno1.g27818.t1"; gene_id "g_27000";

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

    This script is designed for integrating Braker gene preditions onto the gemonic sequence.
    Please use the standard Braker output file!!!

    Usage:
    -g <genome_info>      a text file containing the name and length for each contig/scaffold/chromosome
                          example:    chr1    1000000
                                      chr2    1000000
                                      chr3    1000000
    -i <Braker_file>      the standard Braker output, e.g., braker.gtf
    -r <0-1>              the threshold of overlapping ratio used to merge two mRNA loci. Default: 0.7
    -o <output>           output file
    -h                    print this page

END_OF_TEXT
}
