#!/bin/bash

SEEDS=5
MODELS=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
# orphan_data/ or CASP14/ or CASP13/
for FILE in CASP14/*.fasta;
do
    PREFIX=$(echo ${FILE} | tr -d ".fasta")
    FILENAME=$(echo ${PREFIX} | cut -f2 -d"/")
    echo $PREFIX
    echo $FILENAME
    mkdir -p $PREFIX
    cp $PREFIX.pdb $PREFIX/1org.pdb
    # OmegaFold
    micromamba run -n omegafold omegafold $FILE ${PREFIX}_omega
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference ${PREFIX}_omega/${FILENAME}.pdb --output ${PREFIX}_omega/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # ESMFold
    python fold.py -i $FILE -o ${PREFIX}_esm
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference ${PREFIX}_esm/${FILENAME}.pdb --output ${PREFIX}_esm/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # combine
    mkdir ${PREFIX}_comb
    # copy to comb
    cp ${PREFIX}_esm/${FILENAME}.pdb ${PREFIX}_comb/1esm.pdb
    cp ${PREFIX}_omega/${FILENAME}.pdb ${PREFIX}_comb/1omf.pdb
    # good format for localcolabfold
    mv ${PREFIX}_esm/${FILENAME}.pdb ${PREFIX}_esm/1esm.pdb
    mv ${PREFIX}_omega/${FILENAME}.pdb ${PREFIX}_omega/1omf.pdb
    # Full combined
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}_comb --save-recycles --save-all --num-relax 0 $FILE ${PREFIX}_full_comb
    BEST=$(find ${PREFIX}_full_comb -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_full_comb/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Full ESMFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}_esm --save-recycles --save-all --num-relax 0 $FILE ${PREFIX}_full_esm
    BEST=$(find ${PREFIX}_full_esm -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_full_esm/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Full OmegaFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}_omega --save-recycles --save-all --num-relax 0 $FILE ${PREFIX}_full_omega
    BEST=$(find ${PREFIX}_full_omega -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_full_omega/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Full Truth
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}/ --save-recycles --save-all --overwrite-existing-results --num-relax 0 $FILE ${PREFIX}_full_truth
    BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
    BEST=$(find ${PREFIX}_full_truth -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_full_truth/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Full None
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --num-relax 0 $FILE ${PREFIX}_full_none
    BEST=$(find ${PREFIX}_full_none -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_full_none/scores.json --lddt --local-lddt --tm-score --rigid-scores
    mkdir ${PREFIX}_best/
    cp ${BEST} ${PREFIX}_best/1af1.pdb
    # Single AF
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}_best --save-recycles --save-all --overwrite-existing-results --msa-mode single_sequence --num-relax 0 $FILE ${PREFIX}_single_af
    BEST=$(find ${PREFIX}_single_af -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_single_af/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Single combined
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}_comb --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${PREFIX}_single_comb
    BEST=$(find ${PREFIX}_single_comb -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_single_comb/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Single ESMFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}_esm --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${PREFIX}_single_esm
    BEST=$(find ${PREFIX}_single_esm -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_single_esm/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Single OmegaFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}_omega --save-recycles --save-all --num-relax 0 --msa-mode single_sequence $FILE ${PREFIX}_single_omega
    BEST=$(find ${PREFIX}_single_omega -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_single_omega/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Single Truth
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${PREFIX}/ --save-recycles --save-all --num-relax 0 --msa-mode single_sequence $FILE ${PREFIX}_single_truth
    BEST=$(find ${PREFIX}_single_truth -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_single_truth/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Single None
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --num-relax 0 --msa-mode single_sequence $FILE ${PREFIX}_single_none
    BEST=$(find ${PREFIX}_single_none -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model $PREFIX/1org.pdb --reference $BEST --output ${PREFIX}_single_none/scores.json --lddt --local-lddt --tm-score --rigid-scores
done