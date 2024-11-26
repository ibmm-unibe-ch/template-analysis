#!/bin/bash

function del_dir_and_mk() {
    rm -rf $1
    mkdir -p $1
}
# Taken from https://www.cameo3d.org/modeling/1-month/ in September 2024
CAMEO=("8QJ5A" "8QLUA" "8RIUC" "8RIUD" "8ROHA" "8S4HA" "8VK9D" "8VZID" "8W14A" "8W9UA" "8WC0B" "8WCFA" "8WCGC" "8WDMB" "8WDQA" "8WE5A" "8WEOC" "8WEUA" "8WFHE" "8WKCA" "8WL1A" "8WZCA" "8XYVA" "8YJOA" "8Z4QA" "8ZJ5A" "9AUEA" "9CELD" "9DDYA" "9J91A" "9J9AA")
# Taken from the output of analysis/PDB_hits.ipynb
HITS=("3vsr_A" "7njh_A" "1ytl_A" "3cf4_A" "6zz6_B" "4f9l_D" "6d68_C" "6xn8_A" "2fsx_A" "3l8u_A" "6lhx_A" "6lhy_B" "6z1z_A" "2ogt_A" "6tqo_T" "5yh4_A" "5ocl_A" "3r44_A" "8dce_H" "7v8o_A" "2ew2_A" "3v32_B" "5cwm_A" "1r8w_A" "7y3w_A" "4aty_A" "6cbk_A" "5d3z_A" "4lfe_A" "6ee5_A" "5c2z_A")
OUT_PREFIX=CAMEO1_PDB70
for it in ${!CAMEO[@]}; do
    del_dir_and_mk $OUT_PREFIX/${CAMEO[$it]}
    echo ${CAMEO[$it]}
    echo ${HITS[$it]:0:4}
    echo ${HITS[$it]:5:9999}
    pdb_fetch ${HITS[$it]:0:4} | pdb_selmodel -1 | pdb_selchain -${HITS[$it]:5:9999} | pdb_rplchain -${HITS[$it]:5:9999}:A | pdb_delhetatm | pdb_delinsertion | pdb_tidy | grep ^ATOM >$OUT_PREFIX/${CAMEO[$it]}/1sel.pdb
done
