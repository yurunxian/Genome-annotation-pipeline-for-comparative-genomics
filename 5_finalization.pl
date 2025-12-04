use strict;
use warnings;
use Getopt::Long;

my $gene = "";
my $genome = "";
my $gff_file = "";
my $species = "";
my $output = "";
my $help = "";
my @chr_order = ();
my %gene = ();
my %gene_order = ();
my $usage = "Usage: perl 5_finalization.pl -g <genome_info> -a <gene_gff> --prefix <Species> -o <output>\n";
GetOptions ("g=s" => \$genome,
			"a=s" => \$gff_file,
			"prefix=s" => \$species,
			"o=s" => \$output,
			"h" => \$help);
if ($help) {get_help(); die()}
if (!$gff_file or !$species or !$output) {get_help(); die($usage)}

open (INFO,"<$genome") or die("Please enter the genomic file: -g <genome_info>\n");
open (GENE,"<$gff_file") or die("Please enter the gene gff3 file: -a <gene_gff>\n");
open (OUT,">$output") or die("Please specify the output file: -o <output>\n");

while (<INFO>) {
	chomp;
	my @tmp = split /\s+/,$_;
	push @chr_order,$tmp[0];
}

print "Loading annotated genes... ";
my $gene_number = 0;
while (<GENE>) {
	chomp;
	if (/^#/ or !/\S/) {next}
	my @tmp = split /\t/,$_;
	my $chr = $tmp[0];
	if ($tmp[2] eq "gene" and $tmp[8] =~ /ID=(\S+)/) {
		$gene = $1;
		$gene_number++;
		push @{$gene_order{$chr}}, $gene;
	}
	push @{$gene{$chr}{$gene}}, $_;
}
print " there are $gene_number annotated genes.\n";	

print OUT "##gff-version	3\n";
foreach my $chr (@chr_order) {
	my $number = 0;
	foreach my $gene (@{$gene_order{$chr}}) {
		$number++;
		my $cds = 1;
		my $gene_name = get_number($number,$chr);
		foreach my $line (@{$gene{$chr}{$gene}}) {
			my @tmp = split /\t/, $line;
			if ($tmp[1] eq "transdecoder") {$tmp[1] = "PASA"}
			if ($tmp[2] eq "gene") {$tmp[8] = "ID=$gene_name"}
			elsif ($tmp[2] eq "mRNA") {$tmp[8] = "ID=$gene_name.1;Parent=$gene_name"}
			elsif ($tmp[2] eq "CDS") {
				$tmp[8] = "ID=cds:$gene_name.1:$cds;Parent=$gene_name.1";
				my @exon = @tmp;
				$exon[2] = "exon";
				$exon[8] = "ID=exon:$gene_name.1:$cds;Parent=$gene_name.1";
				print OUT join("\t",@exon),"\n";
				$cds++;
			}
			print OUT join("\t",@tmp),"\n";
		}
		print OUT "###\n";
	}
}

sub get_number {
	my ($number,$chr) = @_;
	my $string = $species;
	if ($chr =~ /(\d+)$/) {
		my $order = $1;
		my $tmp = ("0" x (length(scalar(@chr_order))-length($order))).$order;
		$string .= $tmp;
		$string .= "G";
	}
	my $tmp = ("0" x (4-length($number))).$number;
	$string .= $tmp;
	$string .= "0";
	return $string;
}

sub get_help {
print <<END_OF_TEXT;

    ###########################
    #        Runxian Yu       #
    #  yurunxian\@ibcas.ac.cn  #
    ###########################

    This script is designed for finalizing the gene prediction.
    Please use the output gff3 of 4_integrator.pl!!!

    Usage:
    -g <genome_info>      a text file containing the name and length for each contig/scaffold/chromosome
                          example:    chr1    1000000
                                      chr2    1000000
                                      chr3    1000000
    -a <gene_gff>         the output gff3 files of 4_integrator.pl
    --prefix <species>    the prefix of gene name, e.g., Ath for Arabidopsis thaliana
    -o <output>           output file
    -h                    print this page

END_OF_TEXT
}
