use strict;
use warnings;
use Getopt::Long;

my $gene = "";
my $GeMoMa_file = "";
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
my $usage = "Usage: perl 2_GeMoMa_parser.pl -g <genome_info> -i <GeMoMa_file> -o <output>\n";
GetOptions ("i=s" => \$GeMoMa_file,
			"g=s" => \$genome,
			"r=f" => \$overlap_ratio,
			"o=s" => \$output,
			"h" => \$help) or die(get_help());
if ($help) {get_help(); die()}
if (!$GeMoMa_file or !$genome or !$output) {get_help(); die($usage)}

open (INFO,"<$genome") or die($usage);
open (GEMOMA,"<$GeMoMa_file") or die($usage);
open (OUT,">$output") or die($usage);

print "Loading genome info...";
while (<INFO>) {
	chomp;
	my @tmp = split /\s+/,$_;
	$genome{$tmp[0]} = $tmp[1];
}
print " there are ",scalar keys %genome," contigs.\n";

print "Loading GeMoMa file... ";
my $GeMoMa_number = 0;
while (<GEMOMA>) {
	chomp;
	if (/^#/ or !/\S/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
	if ($tmp[2] eq "gene" and $tmp[8] =~ /ID=(\S+?);/) {
		$gene = $1;
		$number = 1;
		$mRNA = "";
	}
	if ($tmp[2] eq "mRNA" and $tmp[8] =~ /ID=(\S+?);/) {
		$GeMoMa_number++;
		$mRNA = $gene.".".$number;
		push @{$gene{$chr}{$gene}}, $mRNA;
		$strand{$mRNA} = $tmp[6];
		$number++;
	}
	elsif ($tmp[2] eq "CDS" and $tmp[8] =~ /Parent=(\S+)/) {
		$tmp[8] =~ s/Parent=\S+?;/Parent=$mRNA;/g;
		$tmp[8] =~ s/;de=\S+//g;
		$tmp[8] =~ s/;ae=\S+//g;
		my $string = join("\t",@tmp);
		push @{$coding{$mRNA}}, $string;
	}
}
$number = 1;
print " a total of $GeMoMa_number GeMoMa mRNAs are imported!\n";

foreach my $chr (sort keys %genome) {
	print "Processing $chr... ";
	my $anno = "0" x $genome{$chr};
	my $last = "";
	my %anno = ();
	my %gemoma = ();
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
				$gemoma{$mRNA} = $begin;
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
								if ($gemoma{$old} ne "fail") {
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
							$gemoma{$old} = "fail";
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
				$gemoma{$mRNA} = $begin;
				$last = $mRNA;
			}
		}
	}
	# release the memory
	%anno = ();
	my @final = grep {$gemoma{$_} ne "fail"} keys %gemoma;
	print "done!\n";

	#prepare output
	print "Writing gff3 file for $chr... ";
	my $output_number = 0;
	foreach my $mRNA (sort {$gemoma{$a} <=> $gemoma{$b}} @final) {
		$output_number++;
		print OUT "$chr\tGeMoMa\tgene\t",$region{$mRNA},"\t.\t",$strand{$mRNA},"\t.\tID=GeMoMa$number\n";
		print OUT "$chr\tGeMoMa\tmRNA\t",$region{$mRNA},"\t.\t",$strand{$mRNA},"\t.\tID=$mRNA;Parent=GeMoMa$number\n";
		print OUT join("\n",@{$coding{$mRNA}}),"\n\n";
		$number++;
	}
	print "a total of $output_number mRNAs are kept!\n";
}
#chr1	GAF	gene	42494	43313	.	+	.	ID=gene_13363;transcripts=1;complete=1;maxTie=0.0;maxEvidence=3;combinedEvidence=3
#chr1	GeMoMa	mRNA	42494	43313	.	+	.	ID=Scurrula_SPA20777_R0;ref-gene=Scurrula_SPA20777.gene;aa=255;score=595;ce=2;rce=1;tae=0;tde=0;tie=0;minSplitReads=0;tpc=1;minCov=78;avgCov=140.1856;pAA=0.6061;iAA=0.4818;nps=0;start=M;stop=*;evidence=3;Parent=gene_13363;sumWeight=3.0;alternative="Ambore
#chr1	GeMoMa	CDS	31744258	31744407	.	+	0	Parent=Amborella_AmTrH1.05G048800.8.v2.1_R2;de=true;pc=1;minCov=206

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

    This script is designed for integrating GeMoMa gene preditions onto the gemonic sequence.
    Please use the standard GeMoMa output file!!!

    Usage:
    -g <genome_info>      a text file containing the name and length for each contig/scaffold/chromosome
                          example:    chr1    1000000
                                      chr2    1000000
                                      chr3    1000000
    -i <GeMoMa_file>      the standard GeMoMa output, e.g., final_annotation.gff3
    -r <0-1>              the threshold of overlapping ratio used to merge two mRNA loci. Default: 0.7
    -o <output>           output file
    -h                    print this page

END_OF_TEXT
}
