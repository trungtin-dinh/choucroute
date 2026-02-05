---
title: "CHOUCROUTE: Unsupervised classification of coded hyperspectral data using an ensemble of univariate statistical tests"
tags:
  - hyperspectral imaging
  - coded aperture
  - CASSI
  - DD-CASSI
  - unsupervised classification
  - statistical tests
  - MATLAB
authors:
  - name: Trung-Tin Dinh
    orcid: 0009-0000-5293-5201
    affiliation: 1, 2
  - name: Hervé Carfantan
    orcid: 0000-0001-7925-9426
    affiliation: 1
  - name: Antoine Monmayrant
    orcid: 0000-0002-6325-4644
    affiliation: 2
  - name: Simon Lacroix
    orcid: 0000-0001-9988-0981
    affiliation: 2
affiliations:
  - name: IRAP, Universite de Toulouse/CNRS/CNES, France
    index: 1
  - name: LAAS-CNRS, France
    index: 2
date: "2026-02-05"
bibliography: paper.bib
---

## Summary

CHOUCROUTE (Coded Hyperspectral Observations Unsupervised Classification Relying On Univariate Tests Ensemble) is a MATLAB implementation of an unsupervised, region-based classification method designed for *coded* hyperspectral acquisitions produced by snapshot systems of the CASSI family, and in particular the Dual-Disperser CASSI (DD-CASSI) architecture [@Gehm2007; @Hemsley2020calib; @Rouxel2023_ddcassi]. Instead of reconstructing a full hyperspectral cube prior to classification, CHOUCROUTE operates directly on the coded measurements, leveraging statistical tests applied to prediction residuals to detect homogeneous regions and build a label map with associated reference spectra.

## Statement of need

Snapshot coded hyperspectral imagers such as DD-CASSI reduce acquisition time and data volume, while preserving spatial structure, by integrating information across selected wavelengths through programmable coded masks [@Gehm2007; @Hemsley2020calib; @Rouxel2023_ddcassi]. However, many downstream analysis pipelines still rely on reconstructing the full hyperspectral cube before performing classification. This reconstruction step can be computationally demanding and its quality depends strongly on the availability of an accurate segmentation into homogeneous regions, especially when the scene includes mixtures or spectral variability [@Hemsley2022_SA].

CHOUCROUTE addresses this gap by providing an unsupervised classification pipeline tailored to coded data: it detects homogeneous regions directly in the coded domain using an ensemble of univariate statistical tests and progressively grows and merges classes under explicit statistical validation.

## Background and problem setting

A coded acquisition can be modeled as a linear operator applied to the hyperspectral scene, followed by additive noise:

\[
\mathbf{d} = \mathbf{H}\mathbf{o} + \mathbf{n}.
\]

In DD-CASSI-like systems, \(\mathbf{H}\) captures the spatio-spectral filtering induced by the coded mask and spectral integration on the detector, and the number of acquisitions \(S\) is typically much smaller than the number of spectral bands \(W\) [@Hemsley2020calib; @Rouxel2023_ddcassi]. Under a separability assumption, a homogeneous region can be represented by a single reference spectrum modulated by a spatial intensity map, enabling fast estimation of a reference spectrum from coded measurements [@Hemsley2022_SA].

CHOUCROUTE combines: (i) reference spectrum estimation from coded data, (ii) prediction of coded measurements within candidate regions, and (iii) statistical validation of the residuals to decide whether the region is homogeneous and compatible with existing classes.

## Method overview

CHOUCROUTE is structured into three main stages:

1. **Homogeneous region detection**  
   Candidate pixel blocks are selected (guided by spatial structure and intensity heuristics). A reference spectrum is estimated from the candidate set using a fast separability-based procedure [@Hemsley2022_SA]. The corresponding coded data are predicted and residuals are computed. Homogeneity is accepted only if residuals satisfy a set of statistical tests (see below).

2. **Region growing**  
   Starting from a validated seed region, neighboring pixels (or blocks) are iteratively considered for inclusion. Each growth step uses the same predict-and-test principle to ensure that adding pixels preserves statistical homogeneity.

3. **Class merging**  
   When a new homogeneous region is detected, CHOUCROUTE tests whether it is statistically compatible with any existing class. Compatibility is evaluated by forming a balanced pooled set of pixels from the two regions/classes, estimating a reference spectrum, predicting coded data, and applying the statistical tests to the residuals. If validated, the regions are merged; otherwise a new class label is created.

This three-stage design yields a label map and a set of class reference spectra, without requiring full cube reconstruction.

## Statistical validation

CHOUCROUTE assumes (locally) additive, centered Gaussian noise on coded measurements, with unknown variance, which is a common approximation in photon-limited imaging when intensities are sufficiently high [@GoureBrun1997]. The homogeneity decision relies on two complementary types of tests applied to residuals:

- a **mean test** to check that residuals are centered (Student t-test) [@Student1908];
- **normality tests** suited to small sample sizes, including Kolmogorov-Smirnov (and Lilliefors variant), Anderson-Darling, and Shapiro-Wilk [@KS1933; @Lilliefors1967; @AndersonDarling1952; @ShapiroWilk1965].

The test ensemble increases robustness across a range of residual behaviors and limited local sample sizes.

## Software description

The repository provides:

- the main MATLAB entry point `choucroute` implementing the detection/growth/merging pipeline;
- helper functions for statistical testing and region handling;
- a self-contained demo generating synthetic coded data with ground-truth labels, allowing users to run the method without external datasets.

Inputs are coded measurements and the associated coding operators, together with a small number of algorithm parameters (e.g., homogeneity threshold and significance level). Outputs include a label image and estimated reference spectra for the discovered classes.

## Demonstration and expected results

The demo illustrates the full workflow:

- generate a simple multi-class label image,
- define smooth, distinct reference spectra (one per class) and per-pixel intensity factors,
- produce coded measurements under a DD-CASSI-inspired linear model,
- run CHOUCROUTE and visualize the estimated label map and class spectra.

The demo is designed for fast execution and reproducibility, and serves as a minimal test for reviewers and users.

## Limitations

CHOUCROUTE is designed for coded acquisitions where a separability-based local model is meaningful. Performance degrades in strong mixture regions or when intra-class spectral variability violates the single-reference-spectrum assumption. As with other region-based approaches, results depend on parameter settings and on the validity of the local noise model.

## Availability

CHOUCROUTE is released as open-source MATLAB code. The repository includes a demo script and documentation to reproduce a complete run from synthetic data.

## References
[1] M. E. Gehm, et al. *Single-shot compressive spectral imaging with dual-disperser architecture*, Optics Express, 15(21) : 1013–14027, Novembre 2007.
[2] L. Drumetz et al. *Blind Hyperspectral Unmixing Using an Extended Linear Mixing Model to Address Spectral Variability*, IEEE Transactions on Image Processing, 25(8) : 3890-3905, Août 2016.
[3] Shapiro, S. S., and Wilk, M. B. *An analysis of variance test for normality (complete samples)*. Biometrika, 52(3–4), 591–611, 1965.
