# DNN — Deep learning for molecular dynamics

Research code for learning structural differences between molecular-dynamics ensembles using a DenseNet-based classifier, followed by saliency analysis of coordinates associated with class predictions. The workflow includes trajectory preprocessing, feature transformation, NN training, and visualization/interpretation.

## Related publication

Plante, Ambrose, Derek M. Shore, Giulia Morra, George Khelashvili, and Harel Weinstein. ‘A Machine Learning Approach for the Discovery of Ligand-Specific Functional Mechanisms of GPCRs’. Molecules 24, no. 11 (2 June 2019): 2097. https://doi.org/10.3390/molecules24112097.

The archive also contains later, system-specific experiment settings; the checked-in defaults should not be assumed to reproduce every analysis in the publication.

## Workflow and code map

| Stage | Files | Purpose |
| --- | --- | --- |
| Protein preparation | `step_0_getprotein/` | Extract and wrap protein trajectories with VMD/Tcl; concatenate trajectories with CatDCD. |
| Augmentation | `step_1_scramble/` | Randomize protein position and orientation so that the NN is forced to consider invariant features (atom-atom distances, angles, etc). |
| Coordinate extraction | `step_2_get_xyz/` | Write coordinates and class labels to HDF5; choose which atoms to include |
| Classification | `step_3A_ML/` | NN training, validation, and testing. |
| Networks | `nets/` | DenseNet-161-derived architectures, including a variant for smaller atom selections. |
| Features | `trajectory_tools/` | Load coordinates, reshape frames into image-like tensors, split data, and normalize inputs. |
| Visualization | `visualize.ipynb`, `step_3B_visualize/` | post analysis notebooks. |
| Sensitivity | `step_4_sensitivity_analysis/` | Guided-backpropagation saliency for selected classes and samples. |
| Alignment | `align.py` | very fast trajectory-alignment algorithm weighting the most stable atoms the highest for the alignment. |

## usage example for step 2

```bash
python step_2_get_xyz/transform.py protein.psf scrambled.dcd labels.dat coordinates.dat 0
```
