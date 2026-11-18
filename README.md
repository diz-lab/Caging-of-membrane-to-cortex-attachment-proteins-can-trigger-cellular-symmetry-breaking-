# Caging of membrane-to-cortex attachment proteins can trigger cellular symmetry breaking

Analysis code for Dar, Tesoro-Moreno, Palaia et al. *Caging of membrane-to-cortex attachment proteins can trigger cellular symmetry breaking.*

Each folder below is one workflow, with its own README, requirements and license. The READMEs give the method, the input data, the settings, the output and how the results were checked.

## Workflows

| Workflow | What it does | Figures | Language |
|---|---|---|---|
| [`dynamic-micropatterning-analysis`](dynamic-micropatterning-analysis) | Spreading and polarity of single cells after release from a micropattern | Fig 7A–D | Jupyter notebook |
| [`frap-diffusion-analysis`](frap-diffusion-analysis) | FRAP recovery curves, fits, half-times and apparent diffusion coefficients | Fig 6B–D, S7A–C | Jupyter notebook |
| [`3d-mask-intensity-quantification`](3d-mask-intensity-quantification) | Fluorescence intensity in 3D cell or cortex masks, and comparison of conditions | S8E, F, H, I | Jupyter notebooks |
| [`smt-first-passage-analysis`](smt-first-passage-analysis) | Single-molecule tracking: trajectory classification, first passage and free exit times, spider plots | Fig 5G,H | Python scripts |
| [`pore-size-analysis`](pore-size-analysis) | Pore size in the actin cortex, by bandpass filter and watershed | Fig S6O | Python script |
| [`strain-nemaitcs-in-vitro-matlab`](strain-nemaitcs-in-vitro-matlab) | Actin networks on supported lipid bilayers: PIV strain, nematic order, temporal autocorrelation, bundling parameter and pair correlation g(r) | Fig 2F-H, 4C,D, S4D-G | MATLAB scripts |

## Code in each workflow

### Dynamic micropatterning

[`Dynamic_Micropatterning_Analysis.ipynb`](dynamic-micropatterning-analysis/Dynamic_Micropatterning_Analysis.ipynb) starts from one mask stack per cell and makes the radial profiles, single-cell plots, heat maps, polarity index, area, circularity and statistics.

### FRAP

[`FRAP_diffusion_analysis.ipynb`](frap-diffusion-analysis/FRAP_diffusion_analysis.ipynb) starts from the Fiji intensity tables and makes the normalised and average recovery curves, single-exponential fits, half-times, mobile fractions, diffusion coefficients and statistics.

### 3D mask intensity

| Step | Notebook | What it does |
|---|---|---|
| 1 | [`1_quantify_intensity_in_3D_masks.ipynb`](3d-mask-intensity-quantification/1_quantify_intensity_in_3D_masks.ipynb) | z-stack and mask → measurement mask and intensity per z slice, for each cell |
| 2 | [`2_compare_conditions.ipynb`](3d-mask-intensity-quantification/2_compare_conditions.ipynb) | Intensities of all cells → normalised values, statistics and plots |

### Single-molecule tracking

Run the steps in this order. See the [workflow README](smt-first-passage-analysis) for the files to prepare between steps.

| Step | Folder | Script | What it does |
|---|---|---|---|
| 1 | [`01_trajectory_classification`](smt-first-passage-analysis/01_trajectory_classification) | [`trajectory_classification.py`](smt-first-passage-analysis/01_trajectory_classification/trajectory_classification.py) | Computes track features and classifies each track |
| 2 | [`02_median_displacement`](smt-first-passage-analysis/02_median_displacement) | [`sort_mobile_tracks.py`](smt-first-passage-analysis/02_median_displacement/sort_mobile_tracks.py) | Keeps the mobile tracks that pass the step-size and duration filters |
| 2 (optional) | [`02_median_displacement`](smt-first-passage-analysis/02_median_displacement) | [`median_displacement.py`](smt-first-passage-analysis/02_median_displacement/median_displacement.py) | Plots the median displacement of one condition folder |
| 3 | [`03_first_passage_time`](smt-first-passage-analysis/03_first_passage_time) | [`first_passage_time.py`](smt-first-passage-analysis/03_first_passage_time/first_passage_time.py) | First passage time, free exit time, lifetime, end-to-end distance and path length |
| 4 | [`04_bimodal_fpt_fet`](smt-first-passage-analysis/04_bimodal_fpt_fet) | [`bimodal_fpt_fet_analysis.py`](smt-first-passage-analysis/04_bimodal_fpt_fet/bimodal_fpt_fet_analysis.py) | The manuscript first passage and free exit time figures, with automatic choice of fit |
| 5 | [`05_spider_plots`](smt-first-passage-analysis/05_spider_plots) | [`spider_plots.py`](smt-first-passage-analysis/05_spider_plots/spider_plots.py) | Spider plots of the trajectories |

### Pore size

[`pore_size_analysis.py`](pore-size-analysis/pore_size_analysis.py) filters each image, segments the pores with a watershed, and measures the size of every pore.

### Strain and nematics in vitro (MATLAB)

Each analysis runs on any TIFF image or stack. Run [`setup_paths.m`](strain-nemaitcs-in-vitro-matlab/setup_paths.m) once per session, then set the input file, folder or list at the top of a `run_*.m` script. Each stack writes its outputs next to it, plus one summary CSV per folder. See the [workflow README](strain-nemaitcs-in-vitro-matlab) for the methods, parameters and pre-processing.

| Folder | Entry point | What it does |
|---|---|---|
| [`00_shared_utils`](strain-nemaitcs-in-vitro-matlab/00_shared_utils) | – | TIFF reading, CSV writing, figure saving and bounded fitting, used by all analyses |
| [`01_PIV_strain`](strain-nemaitcs-in-vitro-matlab/01_PIV_strain) | [`run_PIV_strain.m`](strain-nemaitcs-in-vitro-matlab/01_PIV_strain/run_PIV_strain.m) | PIV between consecutive frames → strain map, strain over time, maximum strain and strain rate |
| [`02_alignment_nematic_order`](strain-nemaitcs-in-vitro-matlab/02_alignment_nematic_order) | [`run_alignment_nematic.m`](strain-nemaitcs-in-vitro-matlab/02_alignment_nematic_order/run_alignment_nematic.m) | Local filament orientation by FFT, and nematic order parameter per frame |
| [`03_temporal_autocorrelation`](strain-nemaitcs-in-vitro-matlab/03_temporal_autocorrelation) | [`run_temporal_autocorr.m`](strain-nemaitcs-in-vitro-matlab/03_temporal_autocorrelation/run_temporal_autocorr.m) | Temporal autocorrelation of the image, fitted with two exponentials to give the decay times |
| [`04_bundling_parameter`](strain-nemaitcs-in-vitro-matlab/04_bundling_parameter) | [`run_bundling_parameter.m`](strain-nemaitcs-in-vitro-matlab/04_bundling_parameter/run_bundling_parameter.m) | Bundling parameter (standard deviation / mean of the thresholded image) for every frame |
| [`05_pair_correlation_gofr`](strain-nemaitcs-in-vitro-matlab/05_pair_correlation_gofr) | [`run_pair_correlation_gofr.m`](strain-nemaitcs-in-vitro-matlab/05_pair_correlation_gofr/run_pair_correlation_gofr.m) | Pair correlation g(r) of a binary image, fitted with an exponential to give the correlation length |

## Getting started

The Python workflows each install on their own. For example:

```bash
git clone https://github.com/diz-lab/Caging-of-membrane-to-cortex-attachment-proteins-can-trigger-cellular-symmetry-breaking-.git
cd Caging-of-membrane-to-cortex-attachment-proteins-can-trigger-cellular-symmetry-breaking-/frap-diffusion-analysis
pip install -r requirements.txt
```

All Python workflows were tested with Python 3.10.

The MATLAB workflow needs MATLAB with the Image Processing Toolbox (for the nematic order analysis). In MATLAB:

```matlab
run('strain-nemaitcs-in-vitro-matlab/setup_paths.m')
edit run_PIV_strain   % set the input and parameters, then Run
```

## Data

The images and processed data are not in this repository. They are deposited in the BioImage Archive: [BioStudies: https://doi.org/10.6019/S-BIAD4315].

## License

MIT. See the `LICENSE` file in each workflow.

## Contact

Lead contact: Alba Diz-Muñoz (diz@embl.de), European Molecular Biology Laboratory, Heidelberg.
Code: Rubén Tesoro-Moreno (dynamic micropatterning, FRAP, 3D mask intensity), Srishti Dar (single-molecule tracking, pore size) and Zachary Sun (strain and nematics in vitro, MATLAB).

Code development was supported by Claude (Anthropic).
