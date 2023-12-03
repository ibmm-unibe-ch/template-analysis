#!/bin/bash

MONOMERSEEDS=3
MONOMERMODELS=5
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
MULTIMERSEEDS=3
MULTIMERMODELS=2
MULTIMERRECYCLES=5

filearray=("CASP14/T1027.fasta" "CASP14/T1026.fasta" "CASP14/T1028.fasta" "CASP14/T1029.fasta")
#for FILE in dips100/*.fasta;
for FILE in ${filearray[@]} ;
do
    if [[ "$FILE" == *"single"* ]]; then
        continue
    fi
    FILENAME=$(basename $FILE .fasta)  
    PARENT_DIR="$(dirname "$FILE")"
    python get_backbone.py --startpdb ${PARENT_DIR}/${FILENAME}.pdb --targetpdb ${PARENT_DIR}/${FILENAME}/backbone.pdb
    mkdir -p ${PARENT_DIR}_comb_single/${FILENAME}
    mkdir -p ${PARENT_DIR}_af_single/${FILENAME}
    mkdir -p ${PARENT_DIR}_af_single_pred/${FILENAME}
    mkdir -p ${PARENT_DIR}_omega_single/${FILENAME}
    mkdir -p ${PARENT_DIR}_esm_single/${FILENAME}
    mkdir -p ${PARENT_DIR}_comb_single/${FILENAME}
    mkdir -p ${PARENT_DIR}_af_single_pred/${FILENAME}
    mkdir -p ${PARENT_DIR}_fullesm_full/${FILENAME}
    mkdir -p ${PARENT_DIR}_fullaf_full/${FILENAME}
    mkdir -p ${PARENT_DIR}_msaaf_full/${FILENAME}
    mkdir -p ${PARENT_DIR}_comb_full/${FILENAME}
    mkdir -p ${PARENT_DIR}_esm_full/${FILENAME}
    mkdir -p ${PARENT_DIR}_omega_full/${FILENAME}
    mkdir -p ${PARENT_DIR}_truth_full/${FILENAME}
    mkdir -p ${PARENT_DIR}_af_full/${FILENAME}
    MONOMERFILE=${PARENT_DIR}/${FILENAME}_single.fasta
    micromamba run -n omegafold omegafold $MONOMERFILE ${PARENT_DIR}_omega_single/${FILENAME}
    python fold.py -i $MONOMERFILE -o ${PARENT_DIR}_esm_single/${FILENAME}
    colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds $MONOMERSEEDS --num-models $MONOMERMODELS --num-relax 0 $MONOMERFILE ${PARENT_DIR}_af_single_pred/${FILENAME}
    MONOMERCOUNT=0
    for MONOMERFILE in ${PARENT_DIR}_omega_single/${FILENAME}/*.pdb;
    do
        cp ${MONOMERFILE} ${PARENT_DIR}_comb_single/${FILENAME}/${MONOMERCOUNT}omf.pdb
        ((MONOMERCOUNT = MONOMERCOUNT + 1))
    done
    MONOMERCOUNT=0
    for MONOMERFILE in ${PARENT_DIR}_esm_single/${FILENAME}/*.pdb;
    do
        cp ${MONOMERFILE} ${PARENT_DIR}_comb_single/${FILENAME}/${MONOMERCOUNT}esm.pdb
        ((MONOMERCOUNT = MONOMERCOUNT + 1))
    done
    for MONOMERFILE in ${PARENT_DIR}_af_single_pred/${FILENAME}/*.a3m;
    do
        MONOMERNAME=$(basename $MONOMERFILE .a3m)
        BEST=$(find ${PARENT_DIR}_af_single_pred/${FILENAME} -name ${MONOMERNAME}${BESTNAME}|tail -1)
        cp ${BEST} ${PARENT_DIR}_af_single/${FILENAME}/${MONOMERNAME}.pdb
    done
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference ${PARENT_DIR}_HDOCK/${FILENAME}/model_1.pdb --output ${PARENT_DIR}_HDOCK/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    # ESMFold
    python fold.py -i $FILE -o ${PARENT_DIR}_fullesm_full/${FILENAME}
    ESMPDB=$(find ${PARENT_DIR}_fullesm_full/${FILENAME} -name "*pdb"|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference $ESMPDB --output ${PARENT_DIR}_fullesm_full/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    python get_backbone.py --startpdb $ESMPDB --targetpdb ${PARENT_DIR}_fullesm_full/${FILENAME}/backbone.pdb
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/backbone.pdb --reference ${PARENT_DIR}_fullesm_full/${FILENAME}/backbone.pdb --output ${PARENT_DIR}_fullesm_full/${FILENAME}/scores_backbone.json --rigid-scores
    # AlphaFold full, no MSA
    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${PARENT_DIR}_fullaf_full/${FILENAME}
    BEST=$(find ${PARENT_DIR}_fullaf_full/${FILENAME} -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference $BEST --output ${PARENT_DIR}_fullaf_full/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    python get_backbone.py --startpdb $BEST --targetpdb ${PARENT_DIR}_fullaf_full/${FILENAME}/backbone.pdb
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/backbone.pdb --reference ${PARENT_DIR}_fullaf_full/${FILENAME}/backbone.pdb --output ${PARENT_DIR}_fullaf_full/${FILENAME}/scores_backbone.json --rigid-scores
    # AlphaFold full, with MSA
    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --save-recycles --save-all --num-relax 0 $FILE ${PARENT_DIR}_msaaf_full/${FILENAME}
    BEST=$(find ${PARENT_DIR}_msaaf_full/${FILENAME} -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference $BEST --output ${PARENT_DIR}_msaaf_full/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    python get_backbone.py --startpdb $BEST --targetpdb ${PARENT_DIR}_msaaf_full/${FILENAME}/backbone.pdb
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/backbone.pdb --reference ${PARENT_DIR}_msaaf_full/${FILENAME}/backbone.pdb --output ${PARENT_DIR}_msaaf_full/${FILENAME}/scores_backbone.json --rigid-scores
    # Comb templates
    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${PARENT_DIR}_comb_single/${FILENAME} --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${PARENT_DIR}_comb_full/${FILENAME}
    BEST=$(find ${PARENT_DIR}_comb_full/${FILENAME} -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference $BEST --output ${PARENT_DIR}_comb_full/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    python get_backbone.py --startpdb $BEST --targetpdb ${PARENT_DIR}_comb_full/${FILENAME}/backbone.pdb
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/backbone.pdb --reference ${PARENT_DIR}_comb_full/${FILENAME}/backbone.pdb --output ${PARENT_DIR}_comb_full/${FILENAME}/scores_backbone.json --rigid-scores
    # ESMFold templates
    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${PARENT_DIR}_esm_single/${FILENAME} --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${PARENT_DIR}_esm_full/${FILENAME}
    BEST=$(find ${PARENT_DIR}_esm_full/${FILENAME} -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference $BEST --output ${PARENT_DIR}_esm_full/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    python get_backbone.py --startpdb $BEST --targetpdb ${PARENT_DIR}_esm_full/${FILENAME}/backbone.pdb
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/backbone.pdb --reference ${PARENT_DIR}_esm_full/${FILENAME}/backbone.pdb --output ${PARENT_DIR}_esm_full/${FILENAME}/scores_backbone.json --rigid-scores
    # OmegaFold templates
    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${PARENT_DIR}_omega_single/${FILENAME} --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${PARENT_DIR}_omega_full/${FILENAME}
    BEST=$(find ${PARENT_DIR}_omega_full/${FILENAME} -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference $BEST --output ${PARENT_DIR}_omega_full/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    python get_backbone.py --startpdb $BEST --targetpdb ${PARENT_DIR}_omega_full/${FILENAME}/backbone.pdb
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/backbone.pdb --reference ${PARENT_DIR}_omega_full/${FILENAME}/backbone.pdb --output ${PARENT_DIR}_omega_full/${FILENAME}/scores_backbone.json --rigid-scores
    # Truth templates
    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${PARENT_DIR}/${FILENAME}_truth --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${PARENT_DIR}_truth_full/${FILENAME}
    BEST=$(find ${PARENT_DIR}_truth_full/${FILENAME} -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference $BEST --output ${PARENT_DIR}_truth_full/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    python get_backbone.py --startpdb $BEST --targetpdb ${PARENT_DIR}_truth_full/${FILENAME}/backbone.pdb
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/backbone.pdb --reference ${PARENT_DIR}_truth_full/${FILENAME}/backbone.pdb --output ${PARENT_DIR}_truth_full/${FILENAME}/scores_backbone.json --rigid-scores
    # AF templates
    colabfold_batch --recycle-early-stop-tolerance 0 --num-recycle $MULTIMERRECYCLES --overwrite-existing-results --random-seed 6217 --num-seeds $MULTIMERSEEDS --num-models $MULTIMERMODELS --templates --custom-template-path ${PARENT_DIR}_af_single/${FILENAME} --save-recycles --save-all --msa-mode single_sequence --num-relax 0 $FILE ${PARENT_DIR}_af_full/${FILENAME}
    BEST=$(find ${PARENT_DIR}_af_full/${FILENAME} -name $BESTNAME|tail -1)
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/${FILENAME}.pdb --reference $BEST --output ${PARENT_DIR}_af_full/${FILENAME}/scores.json --lddt --local-lddt --tm-score --rigid-scores --interface-scores
    python get_backbone.py --startpdb $BEST --targetpdb ${PARENT_DIR}_af_full/${FILENAME}/backbone.pdb
    podman run --rm -v $(pwd):/home openstructure:latest compare-structures --model ${PARENT_DIR}/${FILENAME}/backbone.pdb --reference ${PARENT_DIR}_af_full/${FILENAME}/backbone.pdb --output ${PARENT_DIR}_af_full/${FILENAME}/scores_backbone.json --rigid-scores
done