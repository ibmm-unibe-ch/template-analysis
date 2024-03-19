#!/bin/bash
source paths.sh
SEEDS=1
MODELS=2
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
RECYCLES=3
NOISE_SIZES=( 1 5 10 )

function del_dir_and_mk(){
    rm -rf $1
    mkdir -p $1
}

function create_rfdiffusion(){
    micromamba run -n SE3nv $RFDIFFUSION_PATH $3 inference.output_prefix=$2 inference.input_pdb=$1 inference.num_designs=1 diffuser.partial_T=$4
    TOP_DIR="$(dirname "$2")"
    FILENAME=$(basename $1 .pdb)
    FASPR_OUTPUT_PATH=${TOP_DIR}_FASPR/$FILENAME
    ATTN_OUTPUT_PATH=${TOP_DIR}_ATTN/$FILENAME
    del_dir_and_mk ${FASPR_OUTPUT_PATH}
    del_dir_and_mk ${ATTN_OUTPUT_PATH}
    $FASPR_PATH -i ${2}_0.pdb -o ${FASPR_OUTPUT_PATH}/1fas.pdb -s ${TOP_DIR:0:6}/$FILENAME.rec.fasta
    $MAXIT_PATH -input ${FASPR_OUTPUT_PATH}/1fas.pdb -output ${FASPR_OUTPUT_PATH}/1fas.cif -o 1    
    micromamba run -n attnpacker python attnpack.py --inputpdb ${FASPR_OUTPUT_PATH}/1fas.pdb --outputpath ${ATTN_OUTPUT_PATH}/1att.pdb --numcores 8
    $MAXIT_PATH -input ${ATTN_OUTPUT_PATH}/1att.pdb -output ${ATTN_OUTPUT_PATH}/1att.cif -o 1
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FASPR_OUTPUT_PATH}/1fas.pdb --reference $1 --output ${FASPR_OUTPUT_PATH}/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${ATTN_OUTPUT_PATH}/1att.pdb --reference $1 --output ${ATTN_OUTPUT_PATH}/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
}

function prepare_af2rank() {
    RESNUMS=$(python count_residues.py --inputpdb $1)
    python get_af2rank.py --startpdb $1 --targetpdb $2/1inp.pdb --length $RESNUMS
    cp $2/1inp.pdb $2/test.pdb
    pdb_sort $2/test.pdb | pdb_tidy | pdb_chain -A>$2/1inp.pdb
    $MAXIT_PATH -input $2/1inp.pdb -output $2/1inp.cif -o 1
}

function af2rank_fold_and_score() {
    del_dir_and_mk af2rank_input_13/
    del_dir_and_mk af2rank_alignments_13/1inp
    del_dir_and_mk af2rank_templates_13/1inp
    TOP_DIR="$(dirname "$3")"
    FILENAME=$(basename $3)
    ONE_FOLDER=${TOP_DIR}_1/${FILENAME}
    TWO_FOLDER=${TOP_DIR}_2/${FILENAME}
    del_dir_and_mk $ONE_FOLDER
    del_dir_and_mk $TWO_FOLDER
    touch af2rank_input_13/1inp.fasta
    echo ">1inp">af2rank_input_13/1inp.fasta
    tail -n +2 $2>>af2rank_input_13/1inp.fasta
    cat af2rank_input_13/1inp.fasta
    if [ $5 == "single_sequence" ]
    then
        cp $2 af2rank_alignments_13/1inp/bfd_uniref_hits.a3m
    else
        touch af2rank_alignments_13/1inp/bfd_uniref_hits.a3m
        echo ">1inp">af2rank_alignments_13/1inp/bfd_uniref_hits.a3m
        LENGTH=$(python count_residues.py --inputpdb $4)
        head -c $LENGTH < /dev/zero | tr '\0' '-'>>af2rank_alignments_13/1inp/bfd_uniref_hits.a3m
    fi
    cat af2rank_alignments_13/1inp/bfd_uniref_hits.a3m
    cp ${1}/* af2rank_templates_13/1inp/
    micromamba run -n openfold_env python /home/jgut/GitHub/openfold/run_pretrained_openfold.py af2rank_input_13 af2rank_templates_13/1inp --use_precomputed_alignments af2rank_alignments_13/ --use_custom_template --uniref90_database_path /data/alphafold_database/uniref90/uniref90.fasta --mgnify_database_path /data/alphafold_database/mgnify/mgy_clusters_2022_05.fa --pdb70_database_path /data/alphafold_database/pdb70/pdb70 --uniclust30_database_path /data/alphafold_database/uniref30/Uniref30_2022_02/UniRef30_2022_02 --bfd_database_path /data/alphafold_database/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --jackhmmer_binary_path ~/micromamba/envs/openfold_env/bin/jackhmmer --hhblits_binary_path ~/micromamba/envs/openfold_env/bin/hhblits --hhsearch_binary_path ~/micromamba/envs/openfold_env/bin/hhsearch --kalign_binary_path ~/micromamba/envs/openfold_env/bin/kalign --output_dir ${ONE_FOLDER} --model_device "cuda:0" --config_preset "model_1_ptm" --skip_relaxation
    micromamba run -n openfold_env python /home/jgut/GitHub/openfold/run_pretrained_openfold.py af2rank_input_13 af2rank_templates_13/1inp --use_precomputed_alignments af2rank_alignments_13/ --use_custom_template --uniref90_database_path /data/alphafold_database/uniref90/uniref90.fasta --mgnify_database_path /data/alphafold_database/mgnify/mgy_clusters_2022_05.fa --pdb70_database_path /data/alphafold_database/pdb70/pdb70 --uniclust30_database_path /data/alphafold_database/uniref30/Uniref30_2022_02/UniRef30_2022_02 --bfd_database_path /data/alphafold_database/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --jackhmmer_binary_path ~/micromamba/envs/openfold_env/bin/jackhmmer --hhblits_binary_path ~/micromamba/envs/openfold_env/bin/hhblits --hhsearch_binary_path ~/micromamba/envs/openfold_env/bin/hhsearch --kalign_binary_path ~/micromamba/envs/openfold_env/bin/kalign --output_dir ${TWO_FOLDER} --model_device "cuda:0" --config_preset "model_2_ptm" --skip_relaxation     
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $ONE_FOLDER/predictions/1inp_model_1_ptm_unrelaxed.pdb --reference $4 --output $ONE_FOLDER/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $TWO_FOLDER/predictions/1inp_model_2_ptm_unrelaxed.pdb --reference $4 --output $TWO_FOLDER/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
    rm af2rank_input_13/1inp.fasta
}

function fold_and_score() {
    echo $1
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path $1 --save-recycles --save-all --num-recycle $RECYCLES --msa-mode single_sequence --num-relax 0 $2 $3
    BEST=$(find $3 -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $BEST --reference $4 --output $3/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
}

# CASP13/ or CASP13/
for FILE in CASP13/*.pdb;
do
    PARENT_DIR="$(dirname "$FILE")"
    FILENAME=$(basename $FILE .pdb)
    echo ${FILENAME}
    PREFIX=${FILE:0:6}
    FLAT_FOLDER_1=${PREFIX}_FLAT_1/$FILENAME
    FLAT_1_AF_FOLDER=${PREFIX}_FLAT_1_AF/$FILENAME
    FLAT_1_AF_RANK_FOLDER=${PREFIX}_FLAT_1_RANK/$FILENAME
    FLAT_1_OF_SS_FOLDER=${PREFIX}_FLAT_1_OF_SS/$FILENAME
    FLAT_1_OF_EM_FOLDER=${PREFIX}_FLAT_1_OF_EMPTY/$FILENAME
    FLAT_FOLDER_2=${PREFIX}_FLAT_2/$FILENAME
    FLAT_2_AF_FOLDER=${PREFIX}_FLAT_2_AF/$FILENAME
    FLAT_2_AF_RANK_FOLDER=${PREFIX}_FLAT_2_RANK/$FILENAME
    FLAT_2_OF_SS_FOLDER=${PREFIX}_FLAT_2_OF_SS/$FILENAME
    FLAT_2_OF_EM_FOLDER=${PREFIX}_FLAT_2_OF_EMPTY/$FILENAME
    NOISE_FOLDER=${PREFIX}_NOISE/$FILENAME
    NOISE_AF_FOLDER=${PREFIX}_NOISE_AF/$FILENAME
    NOISE_AF_RANK_FOLDER=${PREFIX}_NOISE_RANK/$FILENAME
    NOISE_OF_SS_FOLDER=${PREFIX}_NOISE_OF_SS/$FILENAME
    NOISE_OF_EM_FOLDER=${PREFIX}_NOISE_OF_EMPTY/$FILENAME
    FULL_FOLDER=${PREFIX}_full/$FILENAME
    SINGLE_FOLDER=${PREFIX}_single/$FILENAME
    FASTA_FILE=$PARENT_DIR/$FILENAME.rec.fasta
    del_dir_and_mk ${FLAT_FOLDER_1}
    del_dir_and_mk ${FLAT_1_AF_FOLDER}
    del_dir_and_mk ${FLAT_1_AF_RANK_FOLDER}
    del_dir_and_mk ${FLAT_1_OF_SS_FOLDER}
    del_dir_and_mk ${FLAT_1_OF_EM_FOLDER}
    del_dir_and_mk ${FLAT_FOLDER_2}
    del_dir_and_mk ${FLAT_2_AF_FOLDER}
    del_dir_and_mk ${FLAT_2_AF_RANK_FOLDER}
    del_dir_and_mk ${FLAT_2_OF_SS_FOLDER}
    del_dir_and_mk ${FLAT_2_OF_EM_FOLDER}
    del_dir_and_mk ${NOISE_FOLDER}
    del_dir_and_mk ${NOISE_AF_FOLDER}
    del_dir_and_mk ${NOISE_AF_RANK_FOLDER}
    del_dir_and_mk ${NOISE_OF_SS_FOLDER}
    del_dir_and_mk ${NOISE_OF_EM_FOLDER}
    del_dir_and_mk $FULL_FOLDER 
    del_dir_and_mk $SINGLE_FOLDER
    pdb_tofasta $FILE>$FASTA_FILE
    # Noised with RFDIFFUSION
    CONTIG_MAP=$(python count_residues.py --inputpdb $FILE --rfformat 1)
    for NOISE_SIZE in "${NOISE_SIZES[@]}"
    do
        RF_FOLDER=${PREFIX}_RF_${NOISE_SIZE}/$FILENAME
        FASPR_OUTPUT_PATH=${PREFIX}_RF_${NOISE_SIZE}_FASPR/$FILENAME
        ATTN_OUTPUT_PATH=${PREFIX}_RF_${NOISE_SIZE}_ATTN/$FILENAME
        AF_FASPR_OUTPUT_PATH=${PREFIX}_RF_${NOISE_SIZE}_FASPR_AF/$FILENAME
        AF_ATTN_OUTPUT_PATH=${PREFIX}_RF_${NOISE_SIZE}_ATTN_AF/$FILENAME
        RF_AF_RANK_FOLDER=${PREFIX}_RF_${NOISE_SIZE}_RANK/$FILENAME
        RF_OF_SS_FOLDER=${PREFIX}_RF_${NOISE_SIZE}_OF_SS/$FILENAME
        RF_OF_EM_FOLDER=${PREFIX}_RF_${NOISE_SIZE}_OF_EMPTY/$FILENAME
        del_dir_and_mk ${RF_AF_RANK_FOLDER}
        del_dir_and_mk ${RF_OF_SS_FOLDER}
        del_dir_and_mk ${RF_OF_EM_FOLDER}
        create_rfdiffusion $FILE $RF_FOLDER $CONTIG_MAP $NOISE_SIZE
        fold_and_score $FASPR_OUTPUT_PATH $FASTA_FILE $AF_FASPR_OUTPUT_PATH $FILE
        fold_and_score $ATTN_OUTPUT_PATH $FASTA_FILE $AF_ATTN_OUTPUT_PATH $FILE 
        prepare_af2rank $ATTN_OUTPUT_PATH/1att.pdb $RF_AF_RANK_FOLDER
        af2rank_fold_and_score $RF_AF_RANK_FOLDER $FASTA_FILE ${RF_OF_SS_FOLDER} $FILE "single_sequence"
        af2rank_fold_and_score $RF_AF_RANK_FOLDER $FASTA_FILE ${RF_OF_EM_FOLDER} $FILE "empty"
    done
    # Gaussian noise
    python gaussian_noise_pdb.py --startpdb $FILE --targetpdb ${NOISE_FOLDER}/1noi.pdb --std 1
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${NOISE_FOLDER}/1noi.pdb --reference $FILE --output ${NOISE_FOLDER}/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
    fold_and_score $NOISE_FOLDER $FASTA_FILE $NOISE_AF_FOLDER $FILE 
    prepare_af2rank $NOISE_FOLDER/1noi.pdb $NOISE_AF_RANK_FOLDER
    af2rank_fold_and_score $NOISE_AF_RANK_FOLDER $FASTA_FILE ${NOISE_OF_SS_FOLDER} $FILE "single_sequence"
    af2rank_fold_and_score $NOISE_AF_RANK_FOLDER $FASTA_FILE ${NOISE_OF_EM_FOLDER} $FILE "empty"
    # Flattened 2D
    python flatten_pdb.py --startpdb $FILE --targetpdb ${FLAT_FOLDER_2}/1fla.pdb --pca 2
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FLAT_FOLDER_2}/1fla.pdb --reference $FILE --output ${FLAT_FOLDER_2}/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
    fold_and_score $FLAT_FOLDER_2 $FASTA_FILE $FLAT_2_AF_FOLDER $FILE 
    prepare_af2rank $FLAT_FOLDER_2/1fla.pdb $FLAT_2_AF_RANK_FOLDER
    af2rank_fold_and_score $FLAT_2_AF_RANK_FOLDER $FASTA_FILE ${FLAT_2_OF_SS_FOLDER} $FILE "single_sequence"
    af2rank_fold_and_score $FLAT_2_AF_RANK_FOLDER $FASTA_FILE ${FLAT_2_OF_EM_FOLDER} $FILE "empty"
    # Flattened 1D
    python flatten_pdb.py --startpdb $FILE --targetpdb ${FLAT_FOLDER_1}/1fla.pdb --pca 1
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FLAT_FOLDER_1}/1fla.pdb --reference $FILE --output ${FLAT_FOLDER_1}/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
    fold_and_score $FLAT_FOLDER_1 $FASTA_FILE $FLAT_1_AF_FOLDER $FILE 
    prepare_af2rank $FLAT_FOLDER_1/1fla.pdb $FLAT_1_AF_RANK_FOLDER
    af2rank_fold_and_score $FLAT_1_AF_RANK_FOLDER $FASTA_FILE ${FLAT_1_OF_SS_FOLDER} $FILE "single_sequence"
    af2rank_fold_and_score $FLAT_1_AF_RANK_FOLDER $FASTA_FILE ${FLAT_1_OF_EM_FOLDER} $FILE "empty"
    # Full MSA and no template should be covered by bash-packing.sh
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --num-relax 0 $FASTA_FILE $FULL_FOLDER
    BEST=$(find $FULL_FOLDER -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $BEST --reference $FILE --output ${FULL_FOLDER}/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
    # No template
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FASTA_FILE ${SINGLE_FOLDER}
    BEST=$(find $SINGLE_FOLDER -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $BEST --reference $FILE --output ${SINGLE_FOLDER}/ost_scores.json --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
done