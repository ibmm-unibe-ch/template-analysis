#!/bin/bash

SEEDS=1
MODELS=5
BESTNAME10="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
BESTNAME1="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].r1.pdb"
BESTNAME0="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].r0.pdb"

for FILE in CASP14/*.pdb;
do
    FILENAME=$(basename $FILE .pdb)
    PARENT_DIR="$(dirname "$FILE")"
    AF_FULL_FOLDER=${FILE:0:6}_AF_FULL/$FILENAME
    mkdir -p ${AF_FULL_FOLDER}
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --overwrite-existing-results --num-recycle 10 --num-relax 0 $PARENT_DIR/$FILENAME.fasta $AF_FULL_FOLDER
    BEST10=$(find $AF_FULL_FOLDER -name $BESTNAME10|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST10 --output $AF_FULL_FOLDER/ost_scores10.json --lddt --local-lddt --tm-score --rigid-scores
    micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST10 --outputpath $AF_FULL_FOLDER/attn_scores10.json
    # iteration 1
    BEST0=$(find $AF_FULL_FOLDER -name $BESTNAME0|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST0 --output $AF_FULL_FOLDER/ost_scores0.json --lddt --local-lddt --tm-score --rigid-scores
    micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST0 --outputpath $AF_FULL_FOLDER/attn_scores0.json
    # iteration 0
    BEST1=$(find $AF_FULL_FOLDER -name $BESTNAME1|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST1 --output $AF_FULL_FOLDER/ost_scores1.json --lddt --local-lddt --tm-score --rigid-scores
    micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST1 --outputpath $AF_FULL_FOLDER/attn_scores1.json

    
    #FASPR_FOLDER=${FILE:0:6}_FASPR/$FILENAME
    #BACKBONE_FOLDER=${FILE:0:6}_BB/$FILENAME
    #FASPR_PATH=${FASPR_FOLDER}/1fas.pdb
    #AF_FOLDER=${FILE:0:6}_AF/$FILENAME
    #BACKBONE_PATH=${BACKBONE_FOLDER}/1bb1.pdb
    #BB_AF_FOLDER=${FILE:0:6}_BB_AF/$FILENAME
    #mkdir -p ${FASPR_FOLDER}
    #mkdir -p ${AF_FOLDER}
    #mkdir -p ${BACKBONE_FOLDER}
    #mkdir -p ${BB_AF_FOLDER}
    #echo ${FILE}
    #echo ${FILENAME}
    #FASPR/FASPR -i $FILE -o $FASPR_PATH
    #pdb_chain -A $FILE>temp1
    #cp temp1 $FILE
    #pdbfixer $FILE --output=$FILE --replace-nonstandard --add-residues 
    #pdb_chain -A $FASPR_PATH>temp
    #cp temp $FASPR_PATH
    #pdbfixer $FASPR_PATH --output=$FASPR_PATH --replace-nonstandard --add-residues
    #python get_backbone.py --startpdb $FILE --targetpdb $BACKBONE_PATH
    #podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $FASPR_PATH --output $FASPR_FOLDER/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores
    #micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $FASPR_PATH --outputpath $FASPR_FOLDER/attn_scores.json
    ## FASPR AF
    #colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path $FASPR_FOLDER --save-recycles --save-all --overwrite-existing-results --num-recycle 10 --msa-mode single_sequence --num-relax 0 $PARENT_DIR/$FILENAME.fasta $AF_FOLDER
    ## iteration 10
    #BEST10=$(find $AF_FOLDER -name $BESTNAME10|tail -1)
    #podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST10 --output $AF_FOLDER/ost_scores10.json --lddt --local-lddt --tm-score --rigid-scores
    #micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST10 --outputpath $AF_FOLDER/attn_scores10.json
    ## iteration 1
    #BEST0=$(find $AF_FOLDER -name $BESTNAME0|tail -1)
    #podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST0 --output $AF_FOLDER/ost_scores0.json --lddt --local-lddt --tm-score --rigid-scores
    #micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST0 --outputpath $AF_FOLDER/attn_scores0.json
    ## iteration 0
    #BEST1=$(find $AF_FOLDER -name $BESTNAME1|tail -1)
    #podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST1 --output $AF_FOLDER/ost_scores1.json --lddt --local-lddt --tm-score --rigid-scores
    #micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST1 --outputpath $AF_FOLDER/attn_scores1.json

    #### BB AF --not working
    ###colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path $BACKBONE_FOLDER --save-recycles --save-all --overwrite-existing-results --num-recycle 10 --msa-mode single_sequence --num-relax 0 $PARENT_DIR/$FILENAME.fasta $BB_AF_PATH
    #### iteration 10
    ###BEST10=$(find $BB_AF_PATH -name $BESTNAME10|tail -1)
    ###podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST10 --output $BB_AF_PATH/ost_scores10.json --lddt --local-lddt --tm-score --rigid-scores
    ###micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST10 --outputpath $BB_AF_PATH/attn_scores10.json
    #### iteration 1
    ###BEST0=$(find $BB_AF_PATH -name $BESTNAME0|tail -1)
    ###podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST0 --output $BB_AF_PATH/ost_scores0.json --lddt --local-lddt --tm-score --rigid-scores
    ###micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST0 --outputpath $BB_AF_PATH/attn_scores0.json
    #### iteration 0
    ###BEST1=$(find $BB_AF_PATH -name $BESTNAME1|tail -1)
    ###podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST1 --output $BB_AF_PATH/ost_scores1.json --lddt --local-lddt --tm-score --rigid-scores
    ###micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST1 --outputpath $BB_AF_PATH/attn_scores1.json
done
# may need to adapt generate_topology.tcl to correct folder as well!
#vmd -dispdev text -e generate_topology.tcl
#for FILE in CASP14/*.pdb;
#do
#    FILENAME=$(basename $FILE .pdb)  
#    PARENT_DIR="$(dirname "$FILE")"
#    PSF_FOLDER=${FILE:0:6}_BB_PSF/$FILENAME
#    PSF_PATH=${PSF_FOLDER}/1psf.pdb
#    PSF_AF_FOLDER=${FILE:0:6}_PSF_AF/$FILENAME
#    mkdir -p ${PSF_FOLDER}
#    mkdir -p ${PSF_AF_FOLDER}
#    pdbfixer $PSF_PATH --output=$PSF_PATH --replace-nonstandard --add-residues
#    # PSF
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $PSF_PATH --output $PSF_FOLDER/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores
#    micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $PSF_PATH --outputpath $PSF_FOLDER/attn_scores.json
#    # PSF AF
#    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path $PSF_FOLDER --save-recycles --save-all --overwrite-existing-results --num-recycle 10 --msa-mode single_sequence --num-relax 0 $PARENT_DIR/$FILENAME.fasta $PSF_AF_FOLDER
#    # iteration 10
#    BEST10=$(find $PSF_AF_FOLDER -name $BESTNAME10|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST10 --output $PSF_AF_FOLDER/ost_scores10.json --lddt --local-lddt --tm-score --rigid-scores
#    micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST10 --outputpath $PSF_AF_FOLDER/attn_scores10.json
#    # iteration 1
#    BEST0=$(find $PSF_AF_FOLDER -name $BESTNAME0|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST0 --output $PSF_AF_FOLDER/ost_scores0.json --lddt --local-lddt --tm-score --rigid-scores
#    micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST0 --outputpath $PSF_AF_FOLDER/attn_scores0.json
#    # iteration 0
#    BEST1=$(find $PSF_AF_FOLDER -name $BESTNAME1|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $FILE --reference $BEST1 --output $PSF_AF_FOLDER/ost_scores1.json --lddt --local-lddt --tm-score --rigid-scores
#    micromamba run -n attnpacker python score_attnpacker.py --predictedpdb $FILE --targetpdb $BEST1 --outputpath $PSF_AF_FOLDER/attn_scores1.json
#done