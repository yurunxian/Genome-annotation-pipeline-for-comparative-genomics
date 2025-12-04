# Genome-annotation-pipeline-for-comparative-genomics
This pipeline is used to integrate gene predictions from PASA-assembler, GeMoMa and Braker3, and filter out TE-related genes using EDTA TE annotation.

## Prepare your input files
### Transcriptome-based gene prediction
1. Build a comprehensive transcript assembly according to https://github.com/PASApipeline/PASApipeline/wiki/PASA_comprehensive_db .
2. Discard any unreliable transcripts indicated in 'compreh_init_build/compreh_init_build.details', and generate a genome-based gene prediction gff3 file https://github.com/PASApipeline/PASApipeline/wiki/PASA_abinitio_training_sets .
3. The final gff3 file should be something looks like **PASA.transcript.fasta.transdecoder.genome.gff3**.

### Homology-based gene prediction
1. Use GeMoMa for homology-based gene prediction https://jstacs.de/index.php/GeMoMa .
2. I **highly recommand** to provide the RNA-seq mapping file of your own species to GeMoMa, i.e., 'r=MAPPED ERE.m=RNA.sorted.bam' to improve the prediction.
3. The final gff3 file should be **final_annotation.gff**.

### *Ab initio* gene prediction
1. Feed both RNA-seq mapping file (tr) and protein sequences from other species to Braker3 to train the model for your species https://github.com/Gaius-Augustus/BRAKER .
2. The final gtf file should be **braker.gtf**.

### TE annotation
The TE annotation is performed by the EDTA pipeline https://github.com/oushujun/EDTA.

## Integrate the gene predictions
You can now use my scripts to integrate gene predictions from PASA-assembler, GeMoMa and Braker3.
1. Use **1_PASA_parser.pl** to remove redundant transcripts in **PASA.transcript.fasta.transdecoder.genome.gff3**.
```
perl 1_PASA_parser.pl -g genome.len -i PASA.transcript.fasta.transdecoder.genome.gff3 -r 0.7 --remove_single_exon -o PASA.clean.gff3

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

```
2. Use **2_GeMoMa_parser.pl** to remove redundant transcripts in **final_annotation.gff**.
```
perl 2_GeMoMa_parser.pl -g genome.len -i final_annotation.gff -r 0.7 -o GeMoMa.gff3
```
3. Use **3_Braker_parser.pl** to remove redundant transcripts in **braker.gtf**.
```
perl 3_Braker_parser.pl -g genome.len -i braker.gtf -r 0.7 -o Braker.gff3.
```
4. Merge these three gff3 files. Priority order: PASA >= GeMoMa > Braker.
```
perl 4_integrator.pl -g genome.len --est PASA.clean.gff3 --homology GeMoMa.gff3 --denovo Braker.gff3 -r 0.7 -o final.gff3
```
5. Re-name the genes, e.g., Ath1G00010 for the first gene on the first contig.
```
perl 5_finalization.pl

    Usage:
    -g <genome_info>      a text file containing the name and length for each contig/scaffold/chromosome
                          example:    chr1    1000000
                                      chr2    1000000
                                      chr3    1000000
    -a <gene_gff>         the output gff3 files of 4_integrator.pl
    --prefix <species>    the prefix of gene name, e.g., Ath for Arabidopsis thaliana
    -o <output>           output file
    -h                    print this page

Usage: perl 5_finalization.pl -g <genome_info> -a <gene_gff> --prefix <Species> -o <output>
```
6. (Optional) Only keep genes with a < 0.1 (or any value you like) overlapping ratio with annotated TEs.
```
perl 6_TE_filter.pl -g genome.len -t xxx.fasta.mod.EDTA.TEanno.gff3 -r 0.1 -a final.gff3 -o final.TE.gff3
```
   
