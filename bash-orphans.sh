#!/bin/bash

SEEDS=5
MODELS=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"

for FILE in orphan_data/*.fasta;
do
    # OmegaFold
    omegafold $FILE ${FILE:0:16}_omega
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference ${FILE:0:16}_omega/${FILE:12:5}pdb --output ${FILE:0:16}_omega/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # ESMFold
    python fold.py -i $FILE -o ${FILE:0:16}_esm
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference ${FILE:0:16}_esm/${FILE:12:5}pdb --output ${FILE:0:16}_esm/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # combine
    mkdir ${FILE:0:16}_comb
    cp ${FILE:0:16}_esm/${FILE:12:5}pdb ${FILE:0:16}_comb/1esm.pdb
    cp ${FILE:0:16}_omega/${FILE:12:5}pdb ${FILE:0:16}_comb/1omf.pdb
    ## Full combined
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}_comb --save-recycles --save-all --num-relax 0 $FILE ${FILE:0:16}_full_comb
    BEST=$(find ${FILE:0:16}_full_comb -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_full_comb/scores.json --lddt --local-lddt --tm-score --rigid-scores
    ## Full ESMFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}_esm --save-recycles --save-all --num-relax 0 $FILE ${FILE:0:16}_full_esm
    BEST=$(find ${FILE:0:16}_full_esm -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_full_esm/scores.json --lddt --local-lddt --tm-score --rigid-scores
    ## Full OmegaFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}_omega --save-recycles --save-all --num-relax 0 $FILE ${FILE:0:16}_full_omega
    BEST=$(find ${FILE:0:16}_full_omega -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_full_omega/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Full Truth
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}/ --save-recycles --save-all --overwrite-existing-results --num-relax 0 $FILE ${FILE:0:16}_full_truth
    BEST=$(find ${FILE:0:16}_full_truth -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_full_truth/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Full None
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --num-relax 0 $FILE ${FILE:0:16}_full_none
    BEST=$(find ${FILE:0:16}_full_none -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_full_none/scores.json --lddt --local-lddt --tm-score --rigid-scores
    mkdir ${FILE:0:16}_best/
    cp ${BEST} ${FILE:0:16}_best/${FILE:12:5}pdb
    ## Single AF
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}_best --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:16}_single_af
    BEST=$(find ${FILE:0:16}_single_af -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_single_af/scores.json --lddt --local-lddt --tm-score --rigid-scores
    ## Single combined
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}_comb --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:16}_single_comb
    BEST=$(find ${FILE:0:16}_single_comb -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_single_comb/scores.json --lddt --local-lddt --tm-score --rigid-scores
    ## Single ESMFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}_esm --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:16}_single_esm
    BEST=$(find ${FILE:0:16}_single_esm -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_single_esm/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Single OmegaFold
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}_omega --save-recycles --save-all --num-relax 0 --msa-mode single_sequence $FILE ${FILE:0:16}_single_omega
    BEST=$(find ${FILE:0:16}_single_omega -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_single_omega/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Single Truth
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path ${FILE:0:16}/ --save-recycles --save-all --num-relax 0 --msa-mode single_sequence $FILE ${FILE:0:16}_single_truth
    BEST=$(find ${FILE:0:16}_single_truth -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_single_truth/scores.json --lddt --local-lddt --tm-score --rigid-scores
    # Single None
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --num-relax 0 --msa-mode single_sequence $FILE ${FILE:0:16}_single_none
    BEST=$(find ${FILE:0:16}_single_none -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:16}/${FILE:12:5}pdb --reference $BEST --output ${FILE:0:16}_single_none/scores.json --lddt --local-lddt --tm-score --rigid-scores
done