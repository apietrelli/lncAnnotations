#! /bin/sh

#set +o posix

get_file_name_from_path() {
        a=$1
        b=${a%.*}
        echo ${b##*/}
}

######################
# Change log
# - 01/12/2016 : Version 1.0
######################
usage() {
	echo -e '\nbrainArray2Flat - Coversion tool from Brain Array package to Flat file\nDESCRIPTION: This tool is an automatic procedure to convert the "probe_tab" and "desc" files to a Flat file for CDF creation. Several filters will be implemented to select probes and probeset to generate a "UNIQUE" signal.\nUSAGE: brainArray2Flat.sh -d <G.VCFs_DIRECTORY> -o <OUT_FILENAME> [ -i <INTERVAL_LIST> ; -d <INT> ] [-h]\n\nNOTE: Please use ABSOLUTE path for every input file !!!\n'
	exit
}

help(){
        echo "
	Program: brainArray2Flat - Coversion tool from Brain Array package to Flat file
        This pipeline convert Brain Array Custom CDF package files (http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp) to a Flat file.
	OUTPUT:
	1. Flat file

	Database and SW version:

		- Software
		GATK version 	= GATK/3.3.0

		- Variant Annotation
		GATK bundle	= 2.8
		dbsnp		= 138
		HapMap		=
		1000G		=

	Usage:   ExoMe_walker3.sh -d <G.VCFs_DIRECTORY> -o <OUT_FILENAME> [ -i <INTERVAL_LIST> ; -d <INT> ] [-h ]

	Flag:
		-d/--dir        PATH	Path were all the g.vcf files are stored
	        -o/--out        PREFIX  Output prefix name
		-i/--interval   BED3	File containing the regions where the SNP call will be performed (chr,start,end)

	Optional Flag:
		-x/--dp		INT	Depth threshold to call LowCov SNP/INDEL [DEFAULT = 15]
					Es. ALL VARIATION < 15 will be tagged as "LowCov"
	Description:
		-h/--help	flag	Help message with all the database and SW version used by this pipeline

######  NOTE: Please use ABSOLUTE path for every input file !!!  ######
"
}

#My Variables
ARRAY_NAME=""
OUT="NULL"
N_PROBE_TH="4"

# Reading arguments
if [ $# -eq 0 ]; then
	echo -e "\nNO ARGUMENTS GIVEN !!!"
	usage
	exit 1
fi

OPTS=`getopt -o p:d:o:n:h -l probe_tab:,desc:,out:,probe_th:,help -- "$@"`
if [ $? != 0 ]
then
        exit 1
fi

eval set -- "$OPTS"

while [ $# -gt 0 ] ; do
	case "$1"
	in
		-p|--probe_tab) PROBE_TAB="$2"; shift;;
    -d|--desc) DESC="$2"; shift;;
		-o|--out) OUT="$2"; shift;;
    -n|--probe_th) N_PROBE_TH="$2"; shift;;

		-h|--help) help
			exit 1;;
		\?) usage
			exit 1;;
		(--) shift; break;;
		(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
		(*) break;;
	esac
	shift
done

# Test for existing files and correct flags

if ! [ -d $PROBE_TAB ]; then
        echo "ERROR: File $PROBE_TAB does not exist !!!"
	usage
	exit 1
fi
if [ -d $DESC ]; then
        echo "ERROR: File $PROBE_TAB does not exist !!!"
	usage
	exit 1
fi
if [[ $N_PROBE_TH != [0-9]* ]]; then
	echo "ERROR: $N_PROBE_TH is not an integer for -n flag !!!"
	usage
	exit 1
fi
if [ $OUT = "NULL" ]; then
        echo "ERROR: OUT sample name was NOT GIVEN !!!"
	usage
	exit 1
fi
######################## NOT WORKING
# if ! [ -e $INTERVAL_LIST ]; then
#         echo "ERROR: File $INTERVAL_LIST does not exist !!!"
# 	usage
# 	exit 1
# else
# 	INTERVAL_LIST_COL=`head -10 $INTERVAL_LIST | awk 'BEGIN{FS="\t";OFS="\t";TEST="PASS"}{if (NF!=3){TEST="FAIL"}}END{print TEST}'`
# 	if [ $INTERVAL_LIST_COL = "FAIL" ]; then
# 		echo "ERROR: File $INTERVAL_LIST does not seem to ba a BED3 !!!"
# 		echo ""
# 		echo "## Example line:"
# 		head -1 $INTERVAL_LIST
# 		echo "##"
# 		usage
# 		exit 1
# 	fi
# fi
# if [ $INTERVAL_LIST = "NULL" ]; then
#         echo "ERROR: Region interval file (BED3) file was NOT GIVEN !!!"
# 	usage
# 	exit 1
# fi
# if [[ $DP_TH != [0-9]* ]]; then
# 	echo "ERROR: $DP_TH is not an integer for -d flag !!!"
# 	usage
# 	exit 1
# fi
################################### END NOT WORKING
### MAIN
### PIPELINE START ####
## 1- Merging
echo "[Pipeline `date +%c`] START"
echo "[Merging Probe-Desc `date +%c`] Merge annotation and sequence information"
join -t "     " <($DESC sort -k1,1 ) <(sort -k1,1 $PROBE_TAB) > "$OUT"_desc_probe_tab.join.tsv

echo "[Merging Probe-Desc `date +%c`] Select probes mapping DIFFERENT ENSG element"
awk 'BEGIN{FS="\t";OFS="\t"}{split($2,anno,"|");print anno[1],$1,$3"_"$4,$0}' "$OUT"_desc_probe_tab.join.tsv | cut -f 1,3 | sort -k2,2 | uniq | cut -f2 | sort | uniq -d > "$OUT"_probe.dup.id

## 2- Probe filtering ####
echo "[Probe Filtering `date +%c`] Filter ENSG overlapping probes"
echo "[Probe Filtering] Add probe_id to probe_tab file"
awk 'BEGIN{FS="\t";OFS="\t"}{print $2"_"$3,$0}' $PROBE_TAB > "$OUT"_probe_tab.probe_id

echo "[Probe Filtering] Put index (1 - if the probe is overlapping different ENSG) on probe_id duplicated"
awk 'BEGIN{FS="\t";OFS="\t"}{print $0,"1"}' "$OUT"_probe.dup.id > "$OUT"_probe.dup.id.index

echo "[Probe Filtering] Join between "$OUT"_probe_tab.probe_id and probe duplicated file index"
join -t " " <(sort -k1,1 "$OUT"_probe.dup.id.index) <(sort -k1,1 "$OUT"_probe_tab.probe_id) > "$OUT"_probe_tab.probe_id.index

echo "[Probe Filtering] Filter the duplicated one (those with index == 1 in the 8th column)"
awk 'BEGIN{FS="\t";OFS="\t"}{if ($8!="1") {print $0} }' "$OUT"_probe_tab.probe_id.index > "$OUT"_probe_tab.probe_flt

## 3- Probeset filtering ####
echo "[Probeset Filtering `date +%c`] Filter probeset below the probe threshold number"
echo "[Probeset Filtering] Count the number of row, representing the probes, for each probeset"
cut -f2 "$OUT"_probe_tab.probe_flt | sort | uniq -c | sed 's/^ *//;s/ / /' > "$OUT"_probe_tab.probe_flt.probeset-count

echo "[Probeset Filtering] Select only those probeset with $N_PROBE_TH or more probes within"
awk -v n_probe_th=$N_PROBE_TH 'BEGIN{FS="\t";OFS="\t"}{if ($1>=$n_probe_th)print}'  "$OUT"_probe_tab.probe_flt.probeset-count | cut -f2> "$OUT"_probe_tab.probe_flt.probeset_flt.probeset_id

echo "[Probeset Filtering] Filter failed probeset identified in probeset filtering"
join -t " " <(awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$0}' "$OUT"_probe_tab.probe_flt | sort -k1,1) <(sort "$OUT"_probe_tab.probe_flt.probeset_flt.probeset_id) > "$OUT".probe_flt.probeset_flt.probe_tab

## 4- Generate Flat file ####
echo "[Generate Flat file `date +%c`] Filter failed probeset identified in probeset filtering"
join -t " " <(awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$0}' "$OUT".probe_flt.probeset_flt.probe_tab | sort -k1,1) <(awk 'BEGIN{FS="\t";OFS="\t"}{split($2,anno,"|");print $1,anno[1]}' $DESC | sort -k1,1) | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$4,$5,$7,$1,$9}' > "$OUT".flat

echo "FILTERING COMPLETE"

if [ ! -s "$OUT".flat ]; then
  echo "ERROR: Something gone WRONG!!!"
  echo "## Check the INPUT files (probe_tab, desc.txt, chip name) and re-start the pipeline ##"
  usage
  exit 1
else
  touch brainArray2Flat.completed
  # Cleaning
fi
echo "[Pipeline `date +%c`] END"
