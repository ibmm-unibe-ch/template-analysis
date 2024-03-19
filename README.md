# [Dissecting AlphaFolds Capabilities with Limited Sequence Information](https://www.biorxiv.org/content/10.1101/2024.03.14.585076v1)

## Abstract
Protein structure prediction, a fundamental challenge in computational biology, aims to predict a protein's 3D structure from its amino acid sequence. This structure is pivotal for elucidating protein functions, interactions, and driving innovations in drug discovery and enzyme engineering. AlphaFold, a powerful deep learning model, has revolutionized this field by leveraging phylogenetic information from multiple sequence alignments (MSAs) to achieve remarkable accuracy in protein structure prediction. However, a key question remains: how well does AlphaFold understand protein structures? This study investigates AlphaFold's capabilities when relying primarily on high-quality template structures, without the additional information provided by MSAs. By designing experiments that probe local and global structural understanding, we aimed to dissect its dependence on specific features and its ability to handle missing information. Our findings revealed AlphaFold's reliance on sterically valid C-&beta; atoms for correctly interpreting structural templates. Additionally, we observed its remarkable ability to recover 3D structures from certain perturbations. Collectively, these results support the hypothesis that AlphaFold has learned an accurate local biophysical energy function. However, this function seems most effective for local interactions. Our work significantly advances understanding of how deep learning models predict protein structures and provides valuable guidance for researchers aiming to overcome limitations in these models.

## Usage
### Installation
1. `make install-environment`
1. also install [DSSP](https://github.com/cmbi/dssp) into `T_A`
1. `make install-attnpacker`
1. install [OpenFold](https://github.com/rostro36/openfold/tree/main/) in it's own environment `openfold_env`
1. install [RFDiffusion](https://github.com/RosettaCommons/RFdiffusion) into it's own environment 
1. install other utilities if not done before:
    - [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold)
    - [FASPR](https://github.com/tommyhuangthu/FASPR)
    - [MAXIT](https://sw-tools.rcsb.org/apps/MAXIT/binary.html)
    - [OpenStructure](https://git.scicore.unibas.ch/schwede/openstructure/)
    - [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)
1. customize `paths.sh` to your paths

### Reproduce side-chain packing
1. perform installation
1. `bash bash-packing.sh` &rightarrow; you have to adjust for CASP13 or CASP14, also in `generate_topology.tcl`
1. get results by running `Packing_results.ipynb`
1. get SASA analysis by running `SASA_analysis.ipynb`
### Reproduce synthetic backbone refinement
1. perform installation
1. `bash bash-synthetic-backbones.sh` &rightarrow; you have to adjust for CASP13 or CASP14
1. get results by running `Backbone_results.ipynb`

## Contact
If there are questions, please file a [GitHub issue](https://github.com/ibmm-unibe-ch/template-analysis/issues) or send an e-mail to thomas.lemmin@unibe.ch and jannik.gut@unibe.ch.