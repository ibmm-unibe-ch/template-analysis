#!/bin/bash

MODELS=5
SEEDS=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
python aggregation_experiments_manipulation.py
mkdir aggregation_test/half_template_manipulated
mkdir aggregation_test/half_template_rotated
cp -rf aggregation_test/low_template_manipulated/2low.pdb aggregation_test/half_template_manipulated/2low.pdb
cp -rf aggregation_test/up_template_manipulated/2upa.pdb aggregation_test/half_template_manipulated/2upa.pdb
cp -rf aggregation_test/low_template_rotated/2low.pdb aggregation_test/half_template_rotated/2low.pdb
cp -rf aggregation_test/up_template_rotated/2upa.pdb aggregation_test/half_template_rotated/2upa.pdb
# halves
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/half_template_manipulated/ --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_half_manipulated
BEST=$(find aggregation_test/2ubq_half_manipulated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_half_manipulated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/half_template_rotated/ --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_half_rotated
BEST=$(find aggregation_test/2ubq_half_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_half_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
# up
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_manipulated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_up_manipulated
BEST=$(find aggregation_test/2ubq_up_manipulated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_up_manipulated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_rotated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_up_rotated
BEST=$(find aggregation_test/2ubq_up_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_up_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
# low
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/low_template_manipulated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_low_manipulated
BEST=$(find aggregation_test/2ubq_low_manipulated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_low_manipulated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/low_template_rotated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_low_rotated
BEST=$(find aggregation_test/2ubq_low_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_low_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
# full
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/full_template --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_full
BEST=$(find aggregation_test/2ubq_full -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_full/scores.json --lddt --local-lddt --tm-score --rigid-scores
# none
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_none
BEST=$(find aggregation_test/2ubq_none -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_none/scores.json --lddt --local-lddt --tm-score --rigid-scores
