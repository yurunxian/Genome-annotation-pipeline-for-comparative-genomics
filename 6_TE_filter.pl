use strict;
use warnings;
use Getopt::Long;

my $gene = "";
my $gff_file = "";
my $TE_file = "";
my $genome = "";
my $output = "";
my $help = "";
my $overlap_ratio = 0.1;
my %gene_order = ();
my %genome = ();
my %gene = ();
my %gene_number = ();
my %TE = ();
my %region = ();
my $usage = "Usage: perl 6_TE_filter.pl -g <genome_info> -t <TE_gff> -a <gene_gff> -o <output>\n";
GetOptions ("g=s" => \$genome,
			"t=s" => \$TE_file,
			"a=s" => \$gff_file,
			"r=f" => \$overlap_ratio,
			"o=s" => \$output,
			"h" => \$help) or die(get_help($usage));
if ($help) {get_help(); die()}
if (!$gff_file or !$genome or !$TE_file or !$output) {get_help(); die($usage)}

open (INFO,"<$genome") or die("Please enter the genomic file: -g <genome_info>\n");
open (GENE,"<$gff_file") or die("Please enter the gene gff3 file: -a <gene_gff>\n");
open (TE,"<$TE_file") or die("Please enter the TE input file: -t <TE_gff>\n");
open (OUT,">$output") or die("Please specify the output file: -o <output>\n");

while (<INFO>) {
	chomp;
	my @tmp = split /\s+/,$_;
	$genome{$tmp[0]} = $tmp[1];
}

print "Loading annotated TEs... ";
my $TE_number = 0;
while (<TE>) {
	chomp;
	if (/^#/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
#	if ($tmp[2] =~ /LTR/ or $tmp[2] =~ /TIR/ or $tmp[2] =~ /helitron/ or $tmp[2] =~ /LINE/ or $tmp[2] =~ /SINE/) {
	if ($tmp[2] =~ /Gypsy/ or $tmp[2] =~ /Copia/) {
		if ($tmp[4]-$tmp[3]+1 < 50000) {
			$TE_number++;
			my $string = $tmp[3]."\t".$tmp[4];
			push @{$TE{$chr}}, $string;
		}
	}
}
print " there are $TE_number annotated TEs.\n";

print "Loading annotated genes... ";
my $gene_number = 0;
while (<GENE>) {
	chomp;
	if (/^#/ or !/\S/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
	if ($tmp[2] eq "gene" and $tmp[8] =~ /ID=(\S+)/) {
		$gene = $1;
		$gene_number{$chr}++;
		$gene_number++;
		my $string = $tmp[3]."\t".$tmp[4];
		$region{$gene} = $string;
		push @{$gene_order{$chr}}, $gene;
	}
	push @{$gene{$chr}{$gene}}, $_;
}
print " there are $gene_number annotated genes.\n";

foreach my $chr (sort keys %genome) {
	my $TE = "0" x $genome{$chr};

	# map the annotated TEs onto the genomic sequence
	foreach my $i (@{$TE{$chr}}) {
		my ($begin,$end) = split /\t/,$i;
		substr($TE,$begin,$end-$begin+1) = "1" x ($end-$begin+1);
	}

	#prepare output
	print "Writing gff3 file for $chr... ";
	my $output_number = 0;

	# map the annotated genes onto the genomic sequence
	foreach my $gene (@{$gene_order{$chr}}) {
		my ($begin,$end) = split /\t/,$region{$gene};

		# calculate the overlapping ratio between the current gene and TE
		my $string = substr($TE,$begin,$end-$begin+1);
		my @cover = ();
		@cover = $string =~ /1/g;

		# if the overlapping region < 100 bp, then keep it
		if (scalar(@cover)/length($string) < 1-$overlap_ratio) {
			$output_number++;
			print OUT join("\n",@{$gene{$chr}{$gene}}),"\n\n";
		}
	}
	print "a total of $output_number/",$gene_number{$chr}-$output_number," genes pass/fail the TE filtering!\n";
}

sub get_help {
print <<END_OF_TEXT;

    ###########################
    #        Runxian Yu       #
    #  yurunxian\@ibcas.ac.cn  #
    ###########################

    This script is designed for filtering out the annotated genes which are overlapped with TEs.
    Please use the output gff3 of EDTA and 4_integrator.pl!!!

    Usage:
    -g <genome_info>      a text file containing the name and length for each contig/scaffold/chromosome
                          example:    chr1    1000000
                                      chr2    1000000
                                      chr3    1000000
    -t <TE_gff>           the standard EDTA output, e.g., genome.fasta.mod.EDTA.TEanno.gff3
    -a <gene_gff>         the output gff3 files of 4_integrator.pl
    -r <0-1>              the largest overlapping ratio between TE(s) and annotated gene allowed. Default: 0.1
    -o <output>           output file
    -h                    print this page

END_OF_TEXT
}
