# lncAnnotations

Aim: Re-annotate probe from Human gene 1.0 chip for lncRNA discovery

## Note from workflow

### generate fasta file from probe.tsv >>> downloaded from Affy  
```
awk 'BEGIN{FS="\t";OFS="\t"}{print ">probe_"$1"_"$2"_"$3"_"$4"_"$5,$6}' probes.gene10st.tsv | tr '\t' '\n' > HuGene-1_0-st-v1.hg19.probe.mod.fa
```
Some probe in probe.tsv are duplicated because “cover” different transcript clusters

### How many unique probe are within the probe.tsv?
```
cut -f1 HuGene-1_0-st-v1.hg19.probe.tab/HuGene-1_0-st-v1.hg19.probe.tab | sort| uniq | wc
804956
```

### How many probes are in the file?

```
grep ">" HuGene-1_0-st-v1.hg19.probe.mod.fa | wc
861493

#### NOT RUN ####
# bowtie [options]* ebwt.library {-1 mate1 -2 mate2 | --12 read | sequenza ] []
# -f                 query input files are (multi-)FASTA .fa/.mfa
# -k                 report up to good alignments per read (default: 1)
# -m                 suppress all alignments if > <int> exist (def: no limit)
# bowtie -f [FASTA] -k 1 -m 1 ebwt[senza estensione] ### valido per Bowtie 1
# ./bowtie -f ./HuGene-1_0-st-v1.hg19.probe.fa/HuGene-1_0-st-v1.hg19.probe.fa -k 1 -m 1 -v 0 ~/Downloads/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set -S HuGene-1_0-st-v1.hg19.probe.mapped.sam
######
```

#### download index bowtie2 from bowtie website

## Uniquely mapped reads only

* Bowtie mapping

```
./bowtie2-2.2.9/bowtie2 -x ~/Downloads/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -p 4 -U HuGene-1_0-st-v1.hg19.probe.mod.fa -f  > HuGene-1_0-st-v1.hg19.probe.mapped.sam

928313 reads; of these:
  928313 (100.00%) were unpaired; of these:
    88008 (9.48%) aligned 0 times
    716191 (77.15%) aligned exactly 1 time
    124114 (13.37%) aligned >1 times
90.52% overall alignment rate
```

* samtools

```
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
```

## Model for Samtools filtering procedure
### (a.k.a. 2). Filtering those reads with multiple alignment reported (XS: present)

```
samtools view -H output.bam > header.sam
samtools view -F 4 output.bam | grep -v "XS:" | cat header.sam - | \
samtools view -b - > unique.bam
rm header.sam
```

## BEDTools download
```
curl http://bedtools.googlecode.com/files/BEDTools-2.26.0.tar.gz > BEDTools.tar.gz
tar -zxvf BEDTools.tar.gz
cd BEDTools
make
sudo cp bin/* /usr/local/bin/
```
### SEE examples
http://bedtools.readthedocs.io/en/latest/content/example-usage.html

### Transform BAM probes into BED
```
bedtools bamtobed -i HuGene-1_0-st-v1.hg19.probe.mapped.unique.bam > HuGene-1_0-st-v1.hg19.probe.mapped.unique.bed
```

# CDF from Brain array

* 29 novembrer 2016

In http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp guys from University of Michigan, have released a custom CDF based on GENCODE 25 (last updated annotation) containing **protein coding** and **non-coding** transcripts.

We downloaded the entire package for the **TRANSCRIPT** version.

1. go to http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
2. Download GENECODET package
3. Save in Dropbox path:lncrna.annotations/hugene10st_Hs_GENECODET_21.0.0

## Probe filter strategy

check awk split sintax at https://www.gnu.org/software/gawk/manual/gawk.html#String-Functions
some examples of join sintax at http://www.theunixschool.com/2012/01/join-command.html
N.B. "<" within the join command open sub-procedures to be run before joining

```
# Merge annotation and sequence information
join -t "     " <(sort -k1,1 hugene10st_Hs_GENECODET_desc.txt) <(sort -k1,1 hugene10st_Hs_GENECODET_probe_tab) > hugene10st_Hs_GENECODET_desc_probe_tab.join.tsv

awk 'BEGIN{FS="\t";OFS="\t"}{split($2,anno,"|");print anno[1],$1,$3"_"$4"_"$6,$0}' hugene10st_Hs_GENECODET_desc_probe_tab.join.tsv | cut -f 1,3 | sort -k2,2 | uniq | cut -f2 | sort | uniq -d > probe.dup.id
```
