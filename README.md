# CHOUCROUTE

CHOUCROUTE (Coded Hyperspectral Observation Unsupervised Classification Relying On Univariate Tests Ensemble) is a MATLAB function for unsupervised classification of coded hyperspectral data acquired by snapshot compressive imaging systems.

Unlike conventional hyperspectral processing pipelines that rely on full hyperspectral cube reconstruction, CHOUCROUTE operates directly on coded measurements. The method is based on region-level statistical homogeneity tests and does not require ground truth labels or prior knowledge of the number of classes.

---

## Features

- Unsupervised classification of coded hyperspectral data
- Iterative processing with detection, growth, and fusion stages
- Statistical decision based on statistical tests
- Fully implemented in MATLAB

---

## Requirements

- MATLAB (tested with R2024a, compatible with R2021a or later)

The following MATLAB toolboxes are required depending on the selected options:

- **Statistics and Machine Learning Toolbox**
  - `ttest`
  - `kstest`
  - `adtest`

- **Optimization Toolbox**
  - `lsqnonneg`

- **Image Processing Toolbox** (optional)
  - `graythresh`  
  Required only if the panchromatic threshold option `t_pan = "auto"` is used.

The Shapiro–Wilk normality test is supported either via the external `swtest` function (MathWorks File Exchange) or through the lightweight internal implementation `swtestlite` included in the toolbox. No external dependency is required for the default configuration.

---

## Installation

Clone the repository and add it to your MATLAB path:

```matlab
git clone https://github.com/trungtin-dinh/choucroute.git
