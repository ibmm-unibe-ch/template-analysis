#!/bin/bash

function del_dir_and_mk() {
    parentdir="$(dirname "$1")"
    echo $parentdir
    rm -rf $parentdir
    mkdir -p $parentdir
}

TARGET_DB_PATH=TargetDB/targetdb
QUERY_DB_PATH=QueryDB/querydb
RESULT_DB_PATH=ResultDB/resultdb
del_dir_and_mk $TARGET_DB_PATH
del_dir_and_mk $QUERY_DB_PATH
del_dir_and_mk $RESULT_DB_PATH

rm queries.fasta
for FILE in CAMEO1/*.rec.fasta; do
    echo $FILE
    grep "S" $FILE >>queries.fasta
done

# Create mmseqs database for the target db
mmseqs createdb /data/lazars/Jannik_targets/pdb.fasta $TARGET_DB_PATH >TargetDB/log
# Create index for the target database
mmseqs createindex $TARGET_DB_PATH tmp >>TargetDB/log
# Create query database
mmseqs createdb queries.fasta $QUERY_DB_PATH >QueryDB/log
# Run the search code
mmseqs search $QUERY_DB_PATH $TARGET_DB_PATH $RESULT_DB_PATH tmp >ResultDB/log
# Convert result to table
mmseqs convertalis $QUERY_DB_PATH $TARGET_DB_PATH $RESULT_DB_PATH PDB_70_CAMEO.tsv --format-output "query,target,pident,qcov" >>ResultDB/log
#https://mmseqs.com/latest/userguide.pdf
#query, target, percentage of identical, coverage of query
