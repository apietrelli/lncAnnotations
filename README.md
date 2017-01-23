# lncAnnotations

Aim: Re-annotate probe from Human gene 1.0 chip for lncRNA discovery

# CDF from Brain array

* 29 novembrer 2016

In http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp guys from University of Michigan, have released a custom CDF based on GENCODE 25 (last updated annotation) containing **protein coding** and **non-coding** transcripts.

We downloaded the entire package for the **TRANSCRIPT** version.

1. Go to http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
2. Download GENECODET package
3. Save in Dropbox path:lncrna.annotations/hugene10st_Hs_GENECODET_21.0.0

## Probe filter strategy and Flat file creation

We applied 2 filters on the probe/probeset downloaded from Brain Array site.

1. Probe mapping multiple ENSGENE element will be filtered out
2. Probeset with less than 4 probes will be filtered out

To apply those filters and to produce a Flat file for the next procedure, we developed a script named brainArray2Flat.sh

### Example brainArray2Flat.sh command-line for Human Gene 1.0 Array
```
./brainArray2Flat.sh -p hugene10st_Hs_GENECODET_21.0.0/hugene10st_Hs_GENECODET_probe_tab -d hugene10st_Hs_GENECODET_21.0.0/hugene10st_Hs_GENECODET_desc.txt -o test
```
### try brainArray2Flat.sh command-line also for Human Gene 2.0 Array
```
./brainArray2Flat.sh -p ./hugene20st_Hs_GENECODET_probe_tab -d ./hugene20st_Hs_GENECODET_desc.txt -o hugene20st_Hs_GENECODET
```

It will produce a **test.flat** file suitable for falt2cdf script

## Flat to CDF procedure

OLD SCRIPT: Follow instructions in aroma-affymetrix to Convert Affymetrix annotation to Flat file (Perl) at: http://www.aroma-project.org/howtos/create_CDF_from_scratch/
then run perl script
```
shell
sed "s/        /,/" hugene10st_Hs_GENECODET_desc.txt > hugene10st_Hs_GENECODET_desc.csv
perl ./../convertProbesetCSV_differentInput.pl hugene10st_Hs_GENECODET_desc.csv hugene10st_Hs_GENECODET.flt.probe_th.probe_tab hugene10st_Hs_GENECODET.flt.probe_th.flat
```

Go to R and run flat2Cdf

1. Load function from file create.custom.CDF.R  

```
#example: flat2Cdf(file="hjay.r1.flat",chipType="hjay",tag="r1,TC")
#file: assumes header...better perhaps to have ... that passes to read.table?; requires header X, Y
#ucol: unit column
#gcol: group column
#col.class: column classes of file (see read.table); NOTE: needs check that right number?
#splitn: parameter that controls the number of initial chunks that are unwrapped (number of characters of unit names used to keep units together for initial chunks)
#rows:
#cols:

```
2. Run function (be care: original flat2Cdf function put "," as tags >>> change to "."  )
```
library(affxparser)
source("/Users/emagene/Dropbox/codes/github/lncAnnotations/flat2Cdf.R")
flat2Cdf("hugene10st_Hs_GENECODET.flt.probe_th.seq.ENSG_Only.flat", chipType="Gene1.0st.lncrna.genes", tag="v21", col.class=c("character","integer","integer","character","character","character"), xynames=c("X","Y"))
```
3. Make CDF package
```
library(makecdfenv)
make.cdf.package("Gene1.0st.lncrna.genes.v21.cdf", compress = FALSE, species="Homo_sapiens", unlink=TRUE)
```
4. install library then load into R
```
# R CMD build --force gene1.0st.lncrna.genes.v21cdf
# R CMD INSTALL gene1.0st.lncrna.genes.v21cdf_1.50.0.tar.gz
library(gene1.0st.lncrna.genes.v21cdf)
```
5. run affy package & good luck

# Creating description file for master table gene 1.0st

```
cut -f 5 hugene10st_Hs_GENECODET.flt.probe_th.seq.ENSG_Only.NOPERL.flat | sort | uniq | sed "/Group_ID/d" > selected.probesets.arraygene1st.cdf
head -1 hugene10st_Hs_GENECODET_desc.txt > header.desc.txt
join -t "" <(sort -k1,1 hugene10st_Hs_GENECODET_desc.txt) <(sort selected.probesets.arraygene1st.cdf) > hugene10st_Hs_ANNOTATION.provv.txt
cat header.desc.txt hugene10st_Hs_ANNOTATION.provv.txt > hugene10st_Hs_ANNOTATION.txt
rm header.desc.txt hugene10st_Hs_ANNOTATION.provv.txt
```

# Creating description file for master table gene 2.0st

```
cut -f 5 hugene20st_Hs_GENECODET.flat.NEW | sort | uniq | sed "/Group_ID/d" > selected.probesets.arraygene2st.cdf
head -1 hugene20st_Hs_GENECODET_desc.txt > header.desc.txt
join -t "" <(sort -k1,1 hugene20st_Hs_GENECODET_desc.txt) <(sort selected.probesets.arraygene2st.cdf) > hugene20st_Hs_ANNOTATION.provv.txt
cat header.desc.txt hugene20st_Hs_ANNOTATION.provv.txt > hugene20st_Hs_ANNOTATION.txt
rm header.desc.txt hugene20st_Hs_ANNOTATION.provv.txt
```

# Creating BED file for probes (mapping to IGV)

1. create tab delimited file with probeID_X_Y  chrom start end

```
awk 'BEGIN{FS="\t";OFS="\t"}{if (NR>1) print $2,$4,$4+25,$5"_"$6}' hugene10st_Hs_GENECODET_mapping.txt | sort | uniq > hugene10st_Hs_GENECODET_mapping.bed
### for gene 2.0st
awk 'BEGIN{FS="\t";OFS="\t"}{if (NR>1) print $2,$4,$4+25,$5"_"$6}' hugene20st_Hs_GENECODET_mapping.txt | sort | uniq > hugene20st_Hs_GENECODET_mapping.bed
```
2. select only those probes included in (merged) flat file

```
cut -f 1 hugene20st_Hs_GENECODET.flat.NEW | sort | uniq > selected.probes
join -t "" <(awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$0}' hugene20st_Hs_GENECODET_mapping.bed | sort -k1,1) selected.probes | cut -f 2-5 | sort | uniq > hugene20st_Hs_GENECODET_mapping.FILTERED.PROBES.bed
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

### create BED file from Genecode.gtf
```
# download gencode.v25.annotation.gtf from genecode website
# remove headers
sed '/^#/ d' gencode.v25.annotation.gtf > gencode.v25.annotation.noheader.gtf
# add "transcript_id" field (mandatory for bed2gtf)
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.v25.annotation.noheader.gtf > gencode.v25.annotation.noheader.wtrid.gtf
gtf2bed < gencode.v25.annotation.noheader.wtrid.gtf > gencode.v25.annotation.bed


```



# CDF ENSG_Only

create CDF file from original Brainarray selecting only those probes:
- whose Entrez GeneID match at least one ENSG (>>> discard secondary annoations)
- included in the flat file generated from the original file

## Brain array 2 Flat (ENSG_Only)

### Download the brainarray ENSG v21 package from the brainarray site:

* Gene 1.0

http://www.affymetrix.com/estore/browse/products.jsp?productId=131453

- Extract the package into Dropbox directory

### Start the pipeline brainArray2Flat.sh

```
./brainArray2Flat.sh -p ~/Dropbox/lncrna.annotations/hugene10st_Hs_ENSG_21.0.0/hugene10st_Hs_ENSG_probe_tab -d ~/Dropbox/lncrna.annotations/hugene10st_Hs_ENSG_21.0.0/hugene10st_Hs_ENSG_desc.txt -o test
```

### Select unique ENSG that have ENTREZ ID in R
```
library(hugene10sthsentrezgcdf)
x = ls(hugene10sthsentrezgcdf)
library(org.Hs.eg.db)
x1 = gsub("_at","",x)
ensglist = sapply(x1,function(y) unlist(mget(y, org.Hs.egENSEMBL, ifnotfound=NA))[1])
entrezg.sel = as.vector(sapply(selected, function(y) strsplit(y,"\\.")[[1]][1]))
temp = ensglist[!is.na(ensglist)]
write.table(temp, file="~/Dropbox/on.going.papers/lncrna.annotations/ENSG.selected.hugene10st.txt", row.names=F,col.names=F,sep="\t", quote=F)
```

###

```
cut -f5 test.flat | sed 's/_at//' | sort | uniq > ENSG.test_flat.id
awk 'BEGIN{FS="\t";OFS="\t"}{split($5,ensg,"_"); print ensg[1],$0}' test.flat | sort -k1,1 | uniq > ENSG_key.test_flat.join

# How many element? (Probe in probesets)
wc ENSG_key.test_flat.join
599107
# How many genes?
cut -f1 ENSG_key.test_flat.join | sort | uniq | wc
21727
```

### Join flat file with selected ENSG from R

```
join -t "     " ENSG_key.test_flat.join <(sort ENSG.selected.hugene10st.txt) > ENSG_key.test_flat

# How many genes after filters?
cut -f1 ENSG_key.test_flat | sort | uniq | wc
17904

```

### Create a new FLAT file

```
cut -f 2- ENSG_key.test_flat | awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1) {print "Probe_ID","X","Y","Probe_Sequence","Group_ID","Unit_ID"; print $1,$2,$3,$4,$5,$5}else{print $1,$2,$3,$4,$5,$5}}' > ENSG_filter.flat
```

## Run flat2CDF
### REMEMBER: modify rows&cols to 1050 if gene10st, to 1600 if gene20st  








# Mapping approach (Deprecated)

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
