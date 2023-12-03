#!/bin/bash

MODELS=5
SEEDS=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
python aggregation_experiments_manipulation.py
mkdir -p aggregation_test/half_template
mkdir -p aggregation_test/half_template_randomized
mkdir -p aggregation_test/half_template_rotated
mkdir -p aggregation_test/half_template_flattened
mkdir -p aggregation_test/half_template_rfdiff
cp -rf aggregation_test/low_template/2low.pdb aggregation_test/half_template/2low.pdb
cp -rf aggregation_test/up_template/2upa.pdb aggregation_test/half_template/2upa.pdb
cp -rf aggregation_test/low_template_randomized/2low.pdb aggregation_test/half_template_randomized/2low.pdb
cp -rf aggregation_test/up_template_randomized/2upa.pdb aggregation_test/half_template_randomized/2upa.pdb
cp -rf aggregation_test/low_template_rotated/2low.pdb aggregation_test/half_template_rotated/2low.pdb
cp -rf aggregation_test/up_template_rotated/2upa.pdb aggregation_test/half_template_rotated/2upa.pdb
cp -rf aggregation_test/low_template_flattened/2low.pdb aggregation_test/half_template_flattened/2low.pdb
cp -rf aggregation_test/up_template_flattened/2upa.pdb aggregation_test/half_template_flattened/2upa.pdb
cp -rf aggregation_test/low_template_rfdiff/2low.pdb aggregation_test/half_template_rfdiff/2low.pdb
cp -rf aggregation_test/up_template_rfdiff/2upa.pdb aggregation_test/half_template_rfdiff/2upa.pdb
# halves
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/half_template_rotated/ --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_half_rotated
BEST=$(find aggregation_test/2ubq_half_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_half_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/half_template_randomized/ --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_half_randomized
BEST=$(find aggregation_test/2ubq_half_randomized -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_half_randomized/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/half_template_rotated/ --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_half_rotated
BEST=$(find aggregation_test/2ubq_half_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_half_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/half_template_flattened/ --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_half_flattened
BEST=$(find aggregation_test/2ubq_half_flattened -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_half_flattened/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/half_template_rfdiff/ --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_half_rfdiff
BEST=$(find aggregation_test/2ubq_half_rfdiff -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_half_rfdiff/scores.json --lddt --local-lddt --tm-score --rigid-scores
# aggregated
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_rotated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/aggregated_rotated
BEST=$(find aggregation_test/aggregated_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/aggregated_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_randomized --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/aggregated_randomized
BEST=$(find aggregation_test/aggregated_randomized -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/aggregated_randomized/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_rotated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/aggregated_rotated
BEST=$(find aggregation_test/aggregated_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/aggregated_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_flattened --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/aggregated_flattened
BEST=$(find aggregation_test/aggregated_flattened -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/aggregated_flattened/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_rfdiff --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/aggregated_rfdiff
BEST=$(find aggregation_test/aggregated_rfdiff -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/aggregated_rfdiff/scores.json --lddt --local-lddt --tm-score --rigid-scores
# up
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_rotated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_up_rotated
BEST=$(find aggregation_test/2ubq_up_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_up_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_randomized --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_up_randomized
BEST=$(find aggregation_test/2ubq_up_randomized -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_up_randomized/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_rotated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_up_rotated
BEST=$(find aggregation_test/2ubq_up_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_up_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_flattened --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_up_flattened
BEST=$(find aggregation_test/2ubq_up_flattened -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_up_flattened/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/up_template_rfdiff --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_up_rfdiff
BEST=$(find aggregation_test/2ubq_up_rfdiff -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_up_rfdiff/scores.json --lddt --local-lddt --tm-score --rigid-scores
# low
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/low_template_rotated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_low_rotated
BEST=$(find aggregation_test/2ubq_low_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_low_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/low_template_randomized --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_low_randomized
BEST=$(find aggregation_test/2ubq_low_randomized -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_low_randomized/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/low_template_rotated --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_low_rotated
BEST=$(find aggregation_test/2ubq_low_rotated -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_low_rotated/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/low_template_flattened --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_low_flattened
BEST=$(find aggregation_test/2ubq_low_flattened -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_low_flattened/scores.json --lddt --local-lddt --tm-score --rigid-scores
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/low_template_rfdiff --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_low_rfdiff
BEST=$(find aggregation_test/2ubq_low_rfdiff -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_low_rfdiff/scores.json --lddt --local-lddt --tm-score --rigid-scores
# full
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --templates --custom-template-path aggregation_test/full_template --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_full
BEST=$(find aggregation_test/2ubq_full -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_full/scores.json --lddt --local-lddt --tm-score --rigid-scores
# none
colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $SEEDS --num-models $MODELS --save-recycles --save-all --msa-mode single_sequence --num-relax 0 aggregation_test/2ubq.fasta aggregation_test/2ubq_none
BEST=$(find aggregation_test/2ubq_none -name $BESTNAME|tail -1)
podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model aggregation_test/full_template/2ubq.pdb --reference $BEST --output aggregation_test/2ubq_none/scores.json --lddt --local-lddt --tm-score --rigid-scores
