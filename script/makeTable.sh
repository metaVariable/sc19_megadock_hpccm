#!/bin/bash -e
CWD=`pwd`

# default params
OUTPUT_DIR=${OUTPUT_DIR:-"./out"}
TABLE_ITEM_MAX=${TABLE_ITEM_MAX:-"100000"}
TSV_SIZE=${TSV_SIZE:-"50"}
INTERACTIVE=${INTERACTIVE:-"1"}
ENABLE_TSV=${ENABLE_TSV:-"0"}

# args
if [ $# -lt 3 ] ; then
  echo "Usage: $0 [RECEPTOR_MATCH] [LIGAND_MATCH] [TABLE_TITLE]"
  echo "e.g.) $0 \*_r.pdb \*_l.pdb test"
  echo "e.g.) INERACTIVE=0 ENABLE_TSV=1 $0 \*_r.pdb \*_l.pdb test_tsv"
  exit
fi

RECEPTOR_MATCH=$1
LIGAND_MATCH=$2
TABLE_TITLE=$3

RESULT_DIR=${RESULT_DIR:-"$OUTPUT_DIR/$TABLE_TITLE/result"}
TABLE_DIR=${TABLE_DIR:-"${OUTPUT_DIR}/$TABLE_TITLE"}

# print params
echo -e " \n\
RECEPTOR_MATCH = ${RECEPTOR_MATCH} \n\
LIGAND_MATCH   = ${LIGAND_MATCH} \n\
TABLE_TITLE    = ${TABLE_TITLE} \n\
TABLE_DIR      = ${TABLE_DIR=} \n\
RESULT_DIR     = ${RESULT_DIR} \n\
TABLE_ITEM_MAX = ${TABLE_ITEM_MAX} \n\
TSV_SIZE       = ${TSV_SIZE} \n\
INTERACTIVE    = ${INTERACTIVE} \n\
ENABLE_TSV     = ${ENABLE_TSV} \n\
"

receptor_targets=`ls $RECEPTOR_MATCH`
ligand_targets=`ls $LIGAND_MATCH`
common_option="-O"
table_file="$TABLE_DIR/$TABLE_TITLE.table"

receptor_num=`echo $receptor_targets | wc -w`
ligand_num=`echo $ligand_targets | wc -w`
num_pairs=$(( $receptor_num * $ligand_num))

if [ $num_pairs -gt $TABLE_ITEM_MAX ] ; then
  echo "Note: # of total docking pairs is over limit (TABLE_ITEM_MAX=$TABLE_ITEM_MAX) so it is trimmed."
  echo
  num_pairs=$TABLE_ITEM_MAX
fi

echo "# of receptormatch       : $receptor_num"
echo "# of ligand match        : $ligand_num"
echo "# of total docking pairs : $num_pairs"
echo

echo "A table file will be generated at:"
echo "  $table_file"
echo
echo "Please check the above parameters are correctly set."
echo
echo "Note: all relative path will be converted into absolute path."
echo 

if [ $INTERACTIVE -ne 0 ]; then
  read -p "> Start create a docking table: ok? (y/N): " yn
  case "$yn" in [yY]*) ;; *) echo -e "\n abort." ; exit ;; esac
  echo
fi

####################################
# generate table
####################################

# initialize
mkdir -p $TABLE_DIR
mkdir -p $RESULT_DIR
TABLE_DIR=$(cd $TABLE_DIR && pwd)
RESULT_DIR=$(cd $RESULT_DIR && pwd)

table_file="$TABLE_DIR/$TABLE_TITLE.table"

rm -f $table_file
touch $table_file

# write header
echo "TITLE=$TABLE_TITLE" >> "$table_file"
echo "PARAM=-R \$1 -L \$2 -o \$3 $common_option" >> "$table_file"

# set counter
progress_bar='####################'
progress_ratio=${#progress_bar}
count_max=$num_pairs
count=0

# write body
for receptor in $receptor_targets ; do

  r_dirpath=$(cd $(dirname $receptor) && pwd)
  r_filename=$(basename $receptor)
  receptor_path="${r_dirpath}/${r_filename}"

  for ligand in $ligand_targets ; do

    count=$(( count + 1 ))
    if [ $count -gt $num_pairs ] ; then break; fi

    l_dirpath=$(cd $(dirname $ligand) && pwd)
    l_filename=$(basename $ligand)
    ligand_path="${l_dirpath}/${l_filename}"

    # set output filepath
    output_path="${RESULT_DIR}/${r_filename%.*}_${l_filename%.*}.out"

    # write a line
    echo -e "${receptor_path}\t${ligand_path}\t${output_path}" >> "$table_file"

    # show progress
    progress=$(( count * progress_ratio / count_max ))
    echo -ne "\r[\t$count / $count_max ] ${progress_bar:0:$progress}"

  done

  if [ $count -gt $num_pairs ] ; then break; fi

done

echo
echo
echo "Table file is generated at:"
echo "  $table_file"

####################################
# generate tsv
####################################

if [ $ENABLE_TSV -eq 0 ] ; then exit; fi

echo 
echo "> create tsv files ..."
echo

TSV_DIR="${TABLE_DIR}/tsv"
rm -rf ${TSV_DIR}
mkdir -p ${TSV_DIR}

# delete header
sed '1,2d' ${table_file} > "${table_file}.tmp"

# split a table into tsv
opts="--numeric-suffixes --lines=${TSV_SIZE} --suffix-length=5"
split ${opts} "${table_file}.tmp" "${TSV_DIR}/tsv."
rm -f "${table_file}.tmp"

# create table-tsv
tsv_targets=$(ls ${TSV_DIR})
num_tsv=`echo $tsv_targets | wc -w`

# reset table
rm -f "$table_file"
touch "$table_file"

# write header
echo "TITLE=$TABLE_TITLE" >> "$table_file"
echo "PARAM=-I \$1 $common_option" >> "$table_file"

# set counter
progress_bar='####################'
progress_ratio=${#progress_bar}
count_max=$num_tsv
count=0

# write body
for tsv in $tsv_targets ; do

  count=$(( count + 1 ))
  # write a line
  echo -e "${TSV_DIR}/${tsv}" >> "$table_file"

  # show progress
  progress=$(( count * progress_ratio / count_max ))
  echo -ne "\r[\t$count / $count_max ] ${progress_bar:0:$progress}"

done

echo
echo
echo "Table-tsv file is generated at:"
echo "  $table_file"