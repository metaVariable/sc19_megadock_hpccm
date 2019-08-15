#!/bin/bash -e
CWD=`pwd`

# default params
INTERACTIVE=${INTERACTIVE:-"1"}
ENABLE_TSV=${ENABLE_TSV:-"1"}
TABLE_ITEM_MAX=${TABLE_ITEM_MAX:-"100000"}
TSV_SIZE=${TSV_SIZE:-"50"}

# args
if [ $# -lt 4 ] ; then
  echo "Usage: $0 [RELATIVE_ROOT] [INPUT_DIR] [RECEPTOR_MATCH] [LIGAND_MATCH] [TABLE_TITLE]"
  echo "e.g.) $0 . data \*_r.pdb \*_l.pdb test"
  echo "e.g.) RUNTIME_RELATIVE_ROOT=/ $0 . data \*_r.pdb \*_l.pdb test"
  echo "e.g.) INERACTIVE=0 TABLE_ITEM_MAX=100 TSV_SIZE=50 $0 . data \*_r.pdb \*_l.pdb test"
  exit
fi

RELATIVE_ROOT=$1
BASE_INPUT_DIR=${2:-"data"}
RECEPTOR_MATCH=$3
LIGAND_MATCH=$4
TABLE_TITLE=$5
RUNTIME_RELATIVE_ROOT=${RUNTIME_RELATIVE_ROOT:-$RELATIVE_ROOT}

BASE_TABLE_DIR="table"
BASE_OUTPUT_DIR="out"

# print params
echo -e " \n\
# PATH
RELATIVE_ROOT         = ${RELATIVE_ROOT} \n\
RUNTIME_RELATIVE_ROOT = ${RUNTIME_RELATIVE_ROOT} \n\
# INPUT \n\
INPUT_DIR             = ${BASE_INPUT_DIR} \n\
RECEPTOR_MATCH        = ${RECEPTOR_MATCH} \n\
LIGAND_MATCH          = ${LIGAND_MATCH} \n\
TABLE_TITLE           = ${TABLE_TITLE} \n\
# OUTPUT \n\
TABLE_DIR             = ${BASE_TABLE_DIR} \n\
OUTPUT_DIR            = ${BASE_OUTPUT_DIR} \n\
# OPTIONS \n\
TABLE_ITEM_MAX        = ${TABLE_ITEM_MAX} \n\
TSV_SIZE              = ${TSV_SIZE} \n\
INTERACTIVE           = ${INTERACTIVE} \n\
ENABLE_TSV            = ${ENABLE_TSV} \n\
"

INPUT_DIR="${RELATIVE_ROOT}/${BASE_INPUT_DIR}"
TABLE_DIR="${RELATIVE_ROOT}/${BASE_TABLE_DIR}/${TABLE_TITLE}"
OUTPUT_DIR="${RELATIVE_ROOT}/${BASE_OUTPUT_DIR}/${TABLE_TITLE}"

receptor_targets=`ls ${INPUT_DIR}/${RECEPTOR_MATCH}`
ligand_targets=`ls ${INPUT_DIR}/${LIGAND_MATCH}`
common_option="-O"
TABLE_FILE="${TABLE_DIR}/${TABLE_TITLE}.table"

receptor_num=`echo ${receptor_targets} | wc -w`
ligand_num=`echo ${ligand_targets} | wc -w`
num_pairs=$(( ${receptor_num} * ${ligand_num}))

if [ ${num_pairs} -gt ${TABLE_ITEM_MAX} ] ; then
  echo "Note: # of total docking pairs is over limit (TABLE_ITEM_MAX=${TABLE_ITEM_MAX}) so it is trimmed."
  echo
  num_pairs=${TABLE_ITEM_MAX}
fi

echo "# of receptormatch       : ${receptor_num}"
echo "# of ligand match        : ${ligand_num}"
echo "# of total docking pairs : ${num_pairs}"
echo

echo "A table file will be generated at:"
echo "  ${TABLE_FILE}"
echo
echo "Please check the above parameters are correctly set."
echo

if [ ${INTERACTIVE} -ne 0 ]; then
  read -p "> Start create a docking table: ok? (y/N): " yn
  case "$yn" in [yY]*) ;; *) echo -e "\n abort." ; exit ;; esac
  echo
fi

if [ ${INTERACTIVE} -ne 0 ]; then
  echo "> The following directories will be removed: "
  echo "  TABLE_DIR  = ${TABLE_DIR}"
  echo "  OUTPUT_DIR = ${OUTPUT_DIR}"
  read -p "> ok? (y/N): " yn
  case "$yn" in [yY]*) ;; *) echo -e "\n abort." ; exit ;; esac
  echo
fi

echo "> Removing directories ..."
rm -rf ${TABLE_DIR}
rm -rf ${OUTPUT_DIR}

####################################
# generate table
####################################

echo "> Start creating a docking table ..."

# initialize
mkdir -p ${TABLE_DIR}
mkdir -p ${OUTPUT_DIR}

# set RUNTIME_RELATIVE_ROOT
RUNTIME_INPUT_DIR="${RUNTIME_RELATIVE_ROOT}/${BASE_INPUT_DIR}"
RUNTIME_TABLE_DIR="${RUNTIME_RELATIVE_ROOT}/${BASE_TABLE_DIR}/${TABLE_TITLE}"
RUNTIME_OUTPUT_DIR="${RUNTIME_RELATIVE_ROOT}/${BASE_OUTPUT_DIR}/${TABLE_TITLE}"
# replace first "//" to "/" if exists
RUNTIME_INPUT_DIR=${RUNTIME_INPUT_DIR/"//"/"/"}
RUNTIME_TABLE_DIR=${RUNTIME_TABLE_DIR/"//"/"/"}
RUNTIME_OUTPUT_DIR=${RUNTIME_OUTPUT_DIR/"//"/"/"}

# write header
echo "TITLE=${TABLE_TITLE}" >> "${TABLE_FILE}"
echo "PARAM=-R \$1 -L \$2 -o \$3 ${common_option}" >> "${TABLE_FILE}"

# set counter
progress_bar='####################'
progress_ratio=${#progress_bar}
count_max=${num_pairs}
count=0

# write body
for receptor in ${receptor_targets} ; do

  r_filename=$(basename $receptor)
  receptor_path="${RUNTIME_INPUT_DIR}/${r_filename}"

  for ligand in ${ligand_targets} ; do

    count=$(( count + 1 ))
    if [ $count -gt ${num_pairs} ] ; then break; fi

    l_filename=$(basename $ligand)
    ligand_path="${RUNTIME_INPUT_DIR}/${l_filename}"

    # set output filepath
    output_path="${RUNTIME_OUTPUT_DIR}/${r_filename%.*}_${l_filename%.*}.out"

    # write a line
    echo -e "${receptor_path}\t${ligand_path}\t${output_path}" >> "${TABLE_FILE}"

    # show progress
    progress=$(( count * progress_ratio / count_max ))
    echo -ne "\r[\t$count / $count_max ] ${progress_bar:0:$progress}"

  done

  if [ $count -gt $num_pairs ] ; then break; fi

done

echo
echo
echo "Table file is generated at:"
echo "  $TABLE_FILE"

####################################
# generate tsv
####################################

if [ $ENABLE_TSV -eq 0 ] ; then exit; fi

echo 
echo "> creating tsv files ..."
echo

TSV_DIR="${TABLE_DIR}/tsv"
RUNTIME_TSV_DIR="${RUNTIME_TABLE_DIR}/tsv"
RUNTIME_TSV_DIR=${RUNTIME_TSV_DIR/"//"/"/"}

rm -rf ${TSV_DIR}
mkdir -p ${TSV_DIR}

# delete header
sed '1,2d' ${TABLE_FILE} > "${TABLE_FILE}.tmp"

# split a table into tsv
opts="--numeric-suffixes --lines=${TSV_SIZE} --suffix-length=5"
split ${opts} "${TABLE_FILE}.tmp" "${TSV_DIR}/tsv."
rm -f "${TABLE_FILE}.tmp"

# create table-tsv
tsv_targets=$(ls ${TSV_DIR})
num_tsv=`echo $tsv_targets | wc -w`

# reset table
rm -f "$TABLE_FILE"
touch "$TABLE_FILE"

# write header
echo "TITLE=$TABLE_TITLE" >> "$TABLE_FILE"
echo "PARAM=-I \$1 $common_option" >> "$TABLE_FILE"

# set counter
progress_bar='####################'
progress_ratio=${#progress_bar}
count_max=$num_tsv
count=0

# write body
for tsv in $tsv_targets ; do

  count=$(( count + 1 ))
  # write a line
  echo -e "${RUNTIME_TSV_DIR}/${tsv}" >> "$TABLE_FILE"

  # show progress
  progress=$(( count * progress_ratio / count_max ))
  echo -ne "\r[\t$count / $count_max ] ${progress_bar:0:$progress}"

done

echo
echo
echo "Table-tsv file is generated at:"
echo "  $TABLE_FILE"