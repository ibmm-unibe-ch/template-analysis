#!/bin/bash
source paths.sh
SEEDS=1
MODELS=2
RECYCLES=0
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
FASPR_EXEC=/data/jgut/template-analysis/FASPR/FASPR
MAXIT_PATH=/data/jgut/template-analysis/maxit-v10.200-prod-src/bin/maxit
RCSBROOT=/data/jgut/template-analysis/maxit-v10.200-prod-src; export RCSBROOT
PATH="$RCSBROOT/bin:"$PATH; export PATH

function del_dir_and_mk() {
    rm -rf $1
    mkdir -p $1
}

function fold_and_score() {
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path $1 --save-recycles --save-all --overwrite-existing-results --num-recycle $RECYCLES --amber --num-relax 2 --relax-max-iterations 0 --use-gpu-relax --msa-mode single_sequence $2 $3
    BEST=$(find $3 -name $BESTNAME | tail -1)
    score $4 $BEST $3
}

function score() {
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $2 --reference $1 --output $3/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
    micromamba run -n attnpacker python utils/score_attnpacker.py --predictedpdb $2 --targetpdb $1 --outputpath $3/attn_scores.json
    python utils/score.py --parentpath $3 --correctpath $1
}

function clean_pdb() {
    pdb_sort $1 | pdb_tidy | pdb_chain -A >temp_packing_13
    cp temp_packing_13 $1
    pdbfixer $1 --output=$1 --replace-nonstandard --add-residues
    $MAXIT_PATH -input $1 -output ${1%???}cif -o 1
}

for FILE in lemmin/*.pdb; do
    FILENAME=$(basename $FILE .pdb)
    echo ${FILENAME}
    PARENT_DIR="$(dirname "$FILE")"
    PREFIX=${FILE:0:6}
    MIN_FASTA=$PARENT_DIR/$FILENAME.rec.fasta
    AF_FULL_FOLDER=${PREFIX}_AF_FULL/$FILENAME
    TEMPLATE_FOLDER=${PREFIX}_TEMPLATE/$FILENAME
    AF_TEMPLATE_FOLDER=${PREFIX}_TEMPLATE_AF/$FILENAME
    FASPR_FOLDER=${PREFIX}_FASPR/$FILENAME
    FASPR_PATH=${FASPR_FOLDER}/1fas.pdb
    FASPR_AF_FOLDER=${PREFIX}_FASPR_AF/$FILENAME
    ATTNPACK_FOLDER=${PREFIX}_ATTNPACK/$FILENAME
    ATTNPACK_PATH=${ATTNPACK_FOLDER}/1pac.pdb
    ATTNPACK_AF_FOLDER=${PREFIX}_ATTNPACK_AF/$FILENAME
    BACKBONE_FOLDER=${PREFIX}_BB/$FILENAME
    BACKBONE_PATH=${BACKBONE_FOLDER}/1bb1.pdb
    BACKBONE_CB_FOLDER=${PREFIX}_BB_CB/$FILENAME
    BACKBONE_CB_PATH=${BACKBONE_CB_FOLDER}/1bbc.pdb
    CB_DUMMY_FOLDER=${PREFIX}_CB_DUMMY/$FILENAME
    CB_DUMMY_PATH=${CB_DUMMY_FOLDER}/1cbd.pdb
    CB_HEURISTIC_FOLDER=${PREFIX}_CB_HEURISTIC/$FILENAME
    CB_HEURISTIC_PATH=${CB_HEURISTIC_FOLDER}/1cbh.pdb
    del_dir_and_mk ${AF_FULL_FOLDER}
    del_dir_and_mk ${TEMPLATE_FOLDER}
    del_dir_and_mk ${AF_TEMPLATE_FOLDER}
    del_dir_and_mk ${ATTNPACK_FOLDER}
    del_dir_and_mk ${ATTNPACK_AF_FOLDER}
    del_dir_and_mk ${FASPR_AF_FOLDER}
    del_dir_and_mk ${FASPR_FOLDER}
    del_dir_and_mk ${BACKBONE_FOLDER}
    del_dir_and_mk ${BACKBONE_CB_FOLDER}
    del_dir_and_mk ${CB_DUMMY_FOLDER}
    del_dir_and_mk ${CB_HEURISTIC_FOLDER}
    clean_pdb $FILE
    pdb_tofasta $FILE >$MIN_FASTA
    # AlphaFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --overwrite-existing-results --amber --num-relax 2 --relax-max-iterations 0 --use-gpu-relax --num-recycle $RECYCLES $MIN_FASTA $AF_FULL_FOLDER
    BEST=$(find $AF_FULL_FOLDER -name $BESTNAME | tail -1)
    score $FILE $BEST $AF_FULL_FOLDER
    # Full template, single sequence
    cp $FILE $TEMPLATE_FOLDER/1cor.pdb
    fold_and_score $TEMPLATE_FOLDER $MIN_FASTA $AF_TEMPLATE_FOLDER $FILE
    # FASPR
    $FASPR_EXEC -i $FILE -o $FASPR_PATH
    clean_pdb $FASPR_PATH
    score $FILE $FASPR_PATH $FASPR_FOLDER
    # FASPR AF
    fold_and_score $FASPR_FOLDER $MIN_FASTA $FASPR_AF_FOLDER $FILE
    # BB AF
    python utils/get_backbone.py --startpdb $FILE --targetpdb $BACKBONE_PATH
    fold_and_score $BACKBONE_FOLDER $MIN_FASTA $BACKBONE_FOLDER $FILE
    # Attnpack
    micromamba run -n attnpacker python utils/attnpack.py --inputpdb $BACKBONE_PATH --outputpath $ATTNPACK_PATH --numcores 8
    score $FILE $ATTNPACK_PATH $ATTNPACK_FOLDER
    # Attnpack AF
    fold_and_score $ATTNPACK_FOLDER $MIN_FASTA $ATTNPACK_AF_FOLDER $FILE
    # CB backbone
    python utils/get_backbone_cb.py --startpdb $FILE --targetpdb $BACKBONE_CB_PATH
    fold_and_score $BACKBONE_CB_FOLDER $MIN_FASTA $BACKBONE_CB_FOLDER $FILE
    # CB dummy
    python utils/place_cb_dummy.py --startpdb $FILE --targetpdb $CB_DUMMY_PATH
    fold_and_score $CB_DUMMY_FOLDER $MIN_FASTA $CB_DUMMY_FOLDER $FILE
    # CB heuristic
    python utils/place_cb_heuristic.py --startpdb $FILE --targetpdb $CB_HEURISTIC_PATH
    clean_pdb $CB_HEURISTIC_PATH
    fold_and_score $CB_HEURISTIC_FOLDER $MIN_FASTA $CB_HEURISTIC_FOLDER $FILE
done
# may need to adapt generate_topology.tcl to correct folder as well!
vmd -dispdev text -e generate_topology.tcl lemmin
for FILE in lemmin/*.pdb; do
    FILENAME=$(basename $FILE .pdb)
    PARENT_DIR="$(dirname "$FILE")"
    echo $FILENAME
    PREFIX=${FILE:0:6}
    PSF_FOLDER=${PREFIX}_BB_PSF/$FILENAME
    PSF_PATH=${PSF_FOLDER}/1psf.pdb
    PSF_AF_FOLDER=${PREFIX}_PSF_AF/$FILENAME
    del_dir_and_mk ${PSF_AF_FOLDER}
    clean_pdb $PSF_PATH
    # PSF
    score $FILE $PSF_PATH $PSF_FOLDER
    # PSF AF
    fold_and_score $PSF_FOLDER $MIN_FASTA $PSF_AF_FOLDER $FILE
done
