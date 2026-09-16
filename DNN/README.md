# DNN — Deep learning for molecular dynamics

Research code for learning structural differences between molecular-dynamics ensembles using a DenseNet-based classifier, followed by saliency analysis of coordinates associated with class predictions. The workflow connects trajectory preprocessing, coordinate-based feature construction, GPU training, and molecular interpretation.

This portfolio copy preserves the research implementation from [weinsteinlab/DNN](https://github.com/weinsteinlab/DNN). Datasets, trained weights, and the original software environment are not included.

## Related publication

Plante, Ambrose, Derek M. Shore, Giulia Morra, George Khelashvili, and Harel Weinstein. ‘A Machine Learning Approach for the Discovery of Ligand-Specific Functional Mechanisms of GPCRs’. Molecules 24, no. 11 (2 June 2019): 2097. https://doi.org/10.3390/molecules24112097.

The archive also contains later, system-specific experiment settings; the checked-in defaults should not be assumed to reproduce every analysis in the publication.

## Workflow and code map

| Stage | Files | Purpose |
| --- | --- | --- |
| Protein preparation | `step_0_getprotein/` | Extract and wrap protein trajectories with VMD/Tcl; concatenate trajectories with CatDCD. |
| Augmentation | `step_1_scramble/` | Randomize protein position and orientation. |
| Coordinate extraction | `step_2_get_xyz/` | Write coordinates and class labels to HDF5 using MDAnalysis; includes a C-alpha selection variant. |
| Classification | `step_3A_ML/` | Train in stages, reload checkpoints, and evaluate predictions. |
| Networks | `nets/` | DenseNet-161-derived architectures, including a variant for smaller atom selections. |
| Features | `trajectory_tools/` | Load coordinates, reshape frames into image-like tensors, split data, and normalize inputs. |
| Visualization | `visualize.ipynb`, `step_3B_visualize/` | Exploratory analysis notebooks. |
| Sensitivity | `step_4_sensitivity_analysis/` | Guided-backpropagation saliency for selected classes and samples. |
| Alignment | `align.py` | C-alpha distance variability used to weight a trajectory alignment. |

## Algorithm

Each selected frame is represented by Cartesian coordinates with three channels (x, y, z). `traj_to_img` chooses a nearly rectangular grid and reshapes the coordinates for convolution; trailing atoms can be discarded when the atom count does not fit the grid. Training scripts split frame indices into training, validation, and test sets, applying the same assignment to each augmented copy of a frame. Normalization statistics are fitted on training data and applied to the other subsets. A DenseNet-derived network predicts class probabilities, with checkpointing and learning-rate decay during staged training. Guided backpropagation relates predictions to input coordinates; saliency represents model sensitivity rather than direct evidence of causal molecular mechanisms.

## Training stages

1. `ml_new_train.py`: initialize a model on a stride-100 subset.
2. `ml_load_train_stride10.py`: continue from a checkpoint on a stride-10 subset.
3. `ml_load_train_fullset.py`: continue on the full selected dataset.
4. `ml_load_train_regularize.py`: continue with the script's regularization settings.
5. `evaluate.py`: load weights and report prediction accuracy.

Review the user-parameter, model-building, and checkpoint sections in each script. The submission scripts are examples for the original SLURM environment.

## Inputs and dependencies

- Python with NumPy, pandas, h5py, scikit-learn, and MDAnalysis.
- Legacy standalone Keras with a TensorFlow backend. APIs such as `K.tensorflow_backend` and `image_dim_ordering()` are used.
- `keras-vis` for saliency; Jupyter and plotting packages for notebooks.
- VMD/Tcl and CatDCD for preprocessing; SLURM and Conda in submission scripts.

No pinned environment was supplied. Exact compatible versions have not been reconstructed or validated; installing the newest Keras/TensorFlow versions is not a verified setup procedure.

Coordinate conversion uses this positional interface (supply your own input paths):

```bash
python step_2_get_xyz/transform.py protein.psf scrambled.dcd labels.dat coordinates.dat 0
```

The `.dat` outputs are HDF5: `dataset` holds coordinates with shape `(frames, atoms, 3)` and `labels` holds class labels. `transform_CA.py` selects C-alpha atoms. Training defaults expect three classes with `OFS`, `OcS`, and `IFS` filenames; adjust them for your experiment.

## Reproduction notes

- Replace absolute lab paths, Python import paths, Conda environments, SLURM settings, and molecular selections with your configuration. Supply trajectories, topology files, and checkpoints.
- `step_1_scramble/submit.sh` references `extract_and_scramble_protein_notsalient.tcl`, which is absent. The included full-protein script is a different variant; choose the intended atom selection explicitly.
- `evaluate.py` contains an unfilled `%s` checkpoint-directory placeholder. Set the actual directory before evaluation.
- `align.py` uses reference frame 3000 for its final fit and inverse-variance weights. Check trajectory length and zero-variance residues before use.
- The original frame-level random split is retained. Nearby MD frames may be correlated; use trajectory- or time-block-held-out validation to assess generalization beyond this split.

## Portfolio cleanup and provenance

Python formatting is standardized without changing executable syntax trees. Original filenames, algorithms, model configurations, training stages, and notebook code cells are retained. Notebook outputs and execution metadata are cleared for readability; rerun with the original inputs to regenerate figures. Generated Python bytecode caches are excluded. Shell and Tcl scripts retain their original behavior and configuration.

Validation covers Python syntax and equivalence to the uploaded source, notebook source preservation, and shell syntax. End-to-end training and numerical reproduction have not been run because the datasets and legacy environment are not supplied.

Source snapshot: `e547946b0488a3a8d045caa7924a31cc3c064a86` from the supplied `DNN-main.zip`. Original authorship and publication credit are retained. No new license is asserted; the supplied archive contains no license file.
