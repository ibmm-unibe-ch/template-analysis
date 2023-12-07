# [Aggregating and optimising AlphaFold2 structure predictions with templates](https://github.com/ibmm-unibe-ch/template-analysis)

In this work-in-progress project we search for the capabilites of [AlphaFold2](https://github.com/deepmind/alphafold) ([localcolabfold](https://github.com/YoshitakaMo/localcolabfold)) on top of the traditional protein folding.
These experiments also serve as datapoints to find out about the learnt concepts of AlphaFold2, like the suggestion by [Roney & Ovchinnikov](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.129.238101) 
For this, we focussed on the oftentimes overlooked template input and single sequence input.
## ToDo
Among other things, these are still things todo:
- test method with an all gaps input instead of single sequence
- explore results of the DIPS100 dataset
- better code quality (reproducibility and readability)
## Experiments
- Monomer refinement of [CASP14](https://predictioncenter.org/casp14/targetlist.cgi) and [OmegaFold](https://github.com/HeliXonProtein/OmegaFold) hard proteins
- Multimerisation with [CASP14](https://predictioncenter.org/casp14/targetlist.cgi) oligomers
- Protein-protein interaction with [DIPS 100](https://github.com/drorlab/DIPS), like [DiffDock-PP](https://github.com/ketatam/DiffDock-PP)
- Recovery of a backbone diffused with [RFDiffusion](https://github.com/RosettaCommons/RFdiffusion) of [O14933](https://alphafold.ebi.ac.uk/entry/O14933)
- Side-chain packing with [CASP14](https://predictioncenter.org/casp14/targetlist.cgi), like [Attnpacker](https://github.com/MattMcPartlon/AttnPacker)
- Template aggregation with a synthetic test case of [O14933](https://alphafold.ebi.ac.uk/entry/O14933)

## Files
The repository needs better clean-up. We will improve reproducibility soon.
- *bash-\** &rightarrow; main script to make experiment
- *Evaluate_\** &rightarrow; evaluation, where you can also find results
- *Makefile* &rightarrow; central piece for reproducibility
- *figures/* &rightarrow; plots of results
- *tables/* &rightarrow; raw results in table form
- *visualisations/* &rightarrow; other visualisations

## Contact
If there are questions, please file a [GitHub issue](https://github.com/ibmm-unibe-ch/template-analysis/issues) or send an e-mail to thomas.lemmin@unibe.ch and jannik.gut@unibe.ch.
