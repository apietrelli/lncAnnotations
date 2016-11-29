## lncAnnotations

## Note from workflow

## generate fasta file from probe.tsv >>> downloaded from Affx  
awk 'BEGIN{FS="\t";OFS="\t"}{print ">probe_"$1"_"$2"_"$3"_"$4"_"$5,$6}' probes.gene10st.tsv | tr '\t' '\n' > HuGene-1_0-st-v1.hg19.probe.mod.fa

## Some probe in probe.tsv are duplicated because “cover” different transcript clusters
## How many unique probe are within the probe.tsv?
cut -f1 HuGene-1_0-st-v1.hg19.probe.tab/HuGene-1_0-st-v1.hg19.probe.tab | sort| uniq | wc
804956

## How many probes are in the file?
grep ">" HuGene-1_0-st-v1.hg19.probe.mod.fa | wc
861493

```
#### NOT RUN ####
# bowtie [options]* ebwt.library {-1 mate1 -2 mate2 | --12 read | sequenza ] []
# -f                 query input files are (multi-)FASTA .fa/.mfa
# -k                 report up to good alignments per read (default: 1)
# -m                 suppress all alignments if > <int> exist (def: no limit)
# bowtie -f [FASTA] -k 1 -m 1 ebwt[senza estensione] ### valido per Bowtie 1
# ./bowtie -f ./HuGene-1_0-st-v1.hg19.probe.fa/HuGene-1_0-st-v1.hg19.probe.fa -k 1 -m 1 -v 0 ~/Downloads/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set -S HuGene-1_0-st-v1.hg19.probe.mapped.sam
######
```
# >>> dowload index bowtie2 from bowtie website


# Uniquely mapped reads only
1. Bowtie mapping
./bowtie2-2.2.9/bowtie2 -x ~/Downloads/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -p 4 -U HuGene-1_0-st-v1.hg19.probe.mod.fa -f  > HuGene-1_0-st-v1.hg19.probe.mapped.sam

928313 reads; of these:
  928313 (100.00%) were unpaired; of these:
    88008 (9.48%) aligned 0 times
    716191 (77.15%) aligned exactly 1 time
    124114 (13.37%) aligned >1 times
90.52% overall alignment rate


2. samtools
samtools view -Sb HuGene-1_0-st-v1.hg19.probe.mapped.sam | samtools view -F 4 -h - | grep -v "XS:" | samtools view -b - > HuGene-1_0-st-v1.hg19.probe.mapped.unique.bam
samtools flagstat HuGene-1_0-st-v1.hg19.probe.mapped.unique.bam

# -F 4 filter all reads with one or more secondary alignments
# extract XS
# flagstat: statistics relative to bam file 

716191 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
716191 + 0 mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

##### MODEL FOR SAMTOOLS FILTERING PROCEDURE
# (a.k.a. 2). Filtering those reads with multiple alignment reported (XS: present)
# samtools view -H output.bam > header.sam
# samtools view -F 4 output.bam | grep -v "XS:" | cat header.sam - | \
# samtools view -b - > unique.bam
# rm header.sam


###
curl http://bedtools.googlecode.com/files/BEDTools-2.26.0.tar.gz > BEDTools.tar.gz
tar -zxvf BEDTools.tar.gz
cd BEDTools
make
sudo cp bin/* /usr/local/bin/

### SEE examples
http://bedtools.readthedocs.io/en/latest/content/example-usage.html
###


bedtools bamtobed -i HuGene-1_0-st-v1.hg19.probe.mapped.unique.bam > HuGene-1_0-st-v1.hg19.probe.mapped.unique.bed
