# Caging of membrane-to-cortex attachment proteins can trigger cellular symmetry breaking

Analysis code for Dar, Tesoro-Moreno, Palaia et al. *Caging of membrane-to-cortex attachment proteins can trigger cellular symmetry breaking.*

Each folder below is one workflow, with its own README, requirements and license. The READMEs give the method, the input data, the settings, the output and how the results were checked.

## Workflows

| Workflow | What it does | Figures | Language |
|---|---|---|---|
| [`dynamic-micropatterning-analysis`](dynamic-micropatterning-analysis) | Spreading and polarity of single cells after release from a micropattern | Fig 7A–D | Jupyter notebook |
| [`frap-diffusion-analysis`](frap-diffusion-analysis) | FRAP recovery curves, fits, half-times and apparent diffusion coefficients | Fig 6B–D, S7A–C | Jupyter notebook |
| [`3d-mask-intensity-quantification`](3d-mask-intensity-quantification) | Fluorescence intensity in 3D cell or cortex masks, and comparison of conditions | S8E, F, H, I | Jupyter notebooks |
| [`smt-first-passage-analysis`](smt-first-passage-analysis) | Single-molecule tracking: trajectory classification, first passage and free exit times, spider plots | [TO ADD] | Python scripts |
| [`pore-size-analysis`](pore-size-analysis) | Pore size in the actin cortex, by bandpass filter and watershed | [TO ADD] | Python script |

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

## Getting started

Each workflow installs on its own. For example:

```bash
git clone https://github.com/diz-lab/Caging-of-membrane-to-cortex-attachment-proteins-can-trigger-cellular-symmetry-breaking-.git
cd Caging-of-membrane-to-cortex-attachment-proteins-can-trigger-cellular-symmetry-breaking-/frap-diffusion-analysis
pip install -r requirements.txt
```

All workflows were tested with Python 3.10.

## Data

The images and processed data are not in this repository. They are deposited in the BioImage Archive: [TO ADD: accession number].

## License

MIT. See the `LICENSE` file in each workflow.

## Contact

Lead contact: Alba Diz-Muñoz (diz@embl.de), European Molecular Biology Laboratory, Heidelberg.
Code: Rubén Tesoro-Moreno (dynamic micropatterning, FRAP, 3D mask intensity) and Srishti Dar (single-molecule tracking, pore size).

Code development was supported by Claude (Anthropic).
