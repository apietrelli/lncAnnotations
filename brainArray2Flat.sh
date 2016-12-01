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
	echo -e '\nbrainArray2Flat - Coversion tool from Brain Array package to Flat file\nDESCRIPTION: This tool is an automatic procedure to convert the "probe_tab" and "desc" files to a Flat file for CDF creation. Several filters will be implemented to select probes and probeset to generate a "UNIQUE" signal.\nUSAGE: brainArray2Flat.sh -p <FILE.probe_tab> -d <FILE.desc.txt> -o <OUT_FILENAME> [ -n <MIN_PROBES_NUMBER> ] [-h]\n'
	exit
}

help(){
        echo "
	Program: brainArray2Flat - Coversion tool from Brain Array package to Flat file
        This pipeline convert Brain Array Custom CDF package files (http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp) to a Flat file.
  INPUT:
  1. Probe_tab file (*_probe_tab from BrainArray CDF)
  2. Probe set Description file (*_desc.txt from BrainArray CDF)
  3. Min number of probes to select a probeset (Optional)

	OUTPUT:
	1. Flat file

	Usage:   brainArray2Flat.sh -p <FILE.probe_tab> -d <FILE.desc.txt> -o <OUT_FILENAME> [ -n <MIN_PROBES_NUMBER> ] [-h]

	Flag:
		-p/--probe_tab  FILE	Probe tab file from Brain Array CDF package
    -d/--desc       FILE  Description file from Brain Array CDF package
	  -o/--out        PREFIX  Output prefix name

	Optional Flag:
		-n/--probe_th		INT	Min number of probes used as threshold to select a probeset [DEFAULT = 4]
					Es. ALL PROBESET <= 4 will be filtered out
	Description:
		-h/--help	      flag	Help message

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

if ! [ -p $PROBE_TAB ]; then
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
