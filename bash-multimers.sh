#!/bin/bash

MONOMERSEEDS=5
MONOMERMODELS=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
MULTIMERSEEDS=3
MULTIMERMODELS=5
MULTIMERRECYCLES=10
for FILE in casp_oligomers/*.fasta;
do
    if [[ "$FILE" == *"single"* ]]; then
        continue
    fi
#    mkdir ${FILE:0:19}_comb_single
#    mkdir ${FILE:0:19}_af_single
#    MONOMERFILE=${FILE:0:19}_single.fasta
#    omegafold $MONOMERFILE ${FILE:0:19}_omega_single
#    python fold.py -i $MONOMERFILE -o ${FILE:0:19}_esm_single
#    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $MONOMERSEEDS --num-models $MONOMERMODELS --num-relax 0 $MONOMERFILE ${FILE:0:19}_af_single_pred
#    MONOMERCOUNT=0
#    for MONOMERFILE in ${FILE:0:19}_omega_single/*.pdb;
#    do
#        cp ${MONOMERFILE} ${FILE:0:19}_comb_single/${MONOMERCOUNT}omf.pdb
#        ((MONOMERCOUNT = MONOMERCOUNT + 1))
#    done
#    MONOMERCOUNT=0
#    for MONOMERFILE in ${FILE:0:19}_esm_single/*.pdb;
#    do
#        cp ${MONOMERFILE} ${FILE:0:19}_comb_single/${MONOMERCOUNT}esm.pdb
#        ((MONOMERCOUNT = MONOMERCOUNT + 1))
#    done
#    for MONOMERFILE in ${FILE:0:19}_af_single_pred/*.a3m;
#    do
#        BEST=$(find ${FILE:0:19}_af_single_pred -name ${MONOMERFILE:35:4}${BESTNAME}|tail -1)
#        cp ${BEST} ${FILE:0:19}_af_single/${MONOMERFILE:35:4}.pdb
#    done
#    # ESMFold
#    python fold.py -i $FILE -o ${FILE:0:19}_fullesm_full
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:19}/${FILE:15:5}pdb --reference ${FILE:0:19}_fullesm_full/${FILE:15:5}pdb --output ${FILE:0:19}_fullesm_full/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
#    # AlphaFold full, no MSA
#    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:19}_fullaf_full
#    BEST=$(find ${FILE:0:19}_fullaf_full -name $BESTNAME|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:19}/${FILE:15:5}pdb --reference $BEST --output ${FILE:0:19}_fullaf_full/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    # AlphaFold full, with MSA
    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --save-recycles --save-all --num-relax 0 $FILE ${FILE:0:19}_msaaf_full
    BEST=$(find ${FILE:0:19}_msaaf_full -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:19}/${FILE:15:5}pdb --reference $BEST --output ${FILE:0:19}_msaaf_full/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    # Comb templates
#    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${FILE:0:19}_comb_single --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:19}_comb_full
#    BEST=$(find ${FILE:0:19}_comb_full -name $BESTNAME|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:19}/${FILE:15:5}pdb --reference $BEST --output ${FILE:0:19}_comb_full/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
#    # ESMFold templates
#    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${FILE:0:19}_esm_single --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:19}_esm_full
#    BEST=$(find ${FILE:0:19}_esm_full -name $BESTNAME|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:19}/${FILE:15:5}pdb --reference $BEST --output ${FILE:0:19}_esm_full/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
#    # OmegaFold templates
#    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${FILE:0:19}_omega_single --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:19}_omega_full
#    BEST=$(find ${FILE:0:19}_omega_full -name $BESTNAME|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:19}/${FILE:15:5}pdb --reference $BEST --output ${FILE:0:19}_omega_full/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
#    # Truth templates
#    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${FILE:0:19}_truth --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:19}_truth_full
#    BEST=$(find ${FILE:0:19}_truth_full -name $BESTNAME|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:19}/${FILE:15:5}pdb --reference $BEST --output ${FILE:0:19}_truth_full/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
#    # AF templates
#    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${FILE:0:19}_af_single --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${FILE:0:19}_af_full
#    BEST=$(find ${FILE:0:19}_af_full -name $BESTNAME|tail -1)
#    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${FILE:0:19}/${FILE:15:5}pdb --reference $BEST --output ${FILE:0:19}_af_full/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
done