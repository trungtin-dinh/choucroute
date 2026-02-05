---
title: "CHOUCROUTE: Unsupervised classification of coded hyperspectral data using univariate statistical tests"
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

CHOUCROUTE (Coded Hyperspectral Observations Unsupervised Classification Relying On Univariate Tests Ensemble) is a MATLAB implementation of an unsupervised, region-based classification method designed for coded hyperspectral acquisitions produced by snapshot systems of the CASSI family, and in particular the Dual-Disperser CASSI (DD-CASSI) architecture [@Gehm2007; @DunlopGray2016; @Hemsley2020calib; @Rouxel2023_ddcassi]. Instead of reconstructing a full hyperspectral cube prior to classification, CHOUCROUTE operates directly on coded measurements by estimating class reference spectra and validating region homogeneity through statistical tests applied to prediction residuals.

## Statement of need

Snapshot coded hyperspectral imagers such as DD-CASSI reduce acquisition time and data volume by multiplexing spatial and spectral information [@Gehm2007; @DunlopGray2016; @Hemsley2020calib; @Rouxel2023_ddcassi]. In many workflows, classification is performed only after reconstructing a hyperspectral cube, which can be computationally demanding and sensitive to modeling errors and noise. Moreover, supervised classification requires labeled ground truth, which is often unavailable or unreliable in real scenes.

CHOUCROUTE addresses these issues by providing an unsupervised classification pipeline tailored to coded data. It produces a label map and estimated reference spectra without requiring full cube reconstruction and without requiring ground truth labels.

## Background and problem setting

A coded acquisition is modeled as a linear operator applied to the hyperspectral scene, followed by additive noise:

$$
\mathbf{d} = \mathbf{H}\mathbf{o} + \mathbf{n}.
$$

In DD-CASSI-like systems, $\mathbf{H}$ captures the spatio-spectral filtering induced by the coded mask and the spectral integration on the detector, and the number of acquisitions $S$ is typically much smaller than the number of spectral bands $W$ [@Hemsley2020calib; @Rouxel2023_ddcassi]. Under a separability assumption, a homogeneous region can be represented by a single reference spectrum modulated by spatial intensity factors, enabling fast estimation of a reference spectrum from coded measurements [@Hemsley2022_SA].

## Method overview

CHOUCROUTE is an iterative region-based pipeline with three stages: detection of statistically homogeneous seed regions, region growing, and class fusion. At each stage, candidate pixels are evaluated by estimating a reference spectrum (via the separability assumption), predicting coded measurements, and computing residuals between observed and predicted coded data. A region is accepted, expanded, or merged only if the residuals remain statistically consistent with the assumed noise model.

## Statistical validation and options

CHOUCROUTE assumes that coded measurement noise can be locally approximated as additive, centered Gaussian noise with unknown variance, an approximation commonly used in optical measurements when photon counts are sufficiently high [@GoureBrun1997]. The residual validation uses a mean test to check that residuals are centered (Student t-test) [@Student1908], and a Gaussianity test to assess consistency with a normal distribution. The default Gaussianity test in CHOUCROUTE is Shapiro-Wilk [@ShapiroWilk1965]. Alternative Gaussianity tests can be selected through an option, including Kolmogorov-Smirnov (with Lilliefors correction when parameters are estimated) [@Kolmogorov1933; @Lilliefors1967] and Anderson-Darling [@AndersonDarling1952].

## Software description

The repository provides the main MATLAB entry point `choucroute`, together with supporting routines for region handling, reference spectrum estimation, prediction, and statistical testing. The algorithm exposes a small set of parameters. In the default configuration, only the intra-class variability threshold $T_\psi$ requires manual tuning, while other parameters are set from data dimensions and standard statistical choices.

## Demonstration and expected results

A self-contained demo script is included to generate a simple synthetic coded dataset with ground-truth regions, smooth class spectra, and per-pixel variability coefficients. The script runs CHOUCROUTE on these coded measurements and produces a label map and estimated class spectra, allowing users and reviewers to verify that the software executes correctly without external data.

## Limitations

CHOUCROUTE is designed for coded acquisitions where a separability-based local model is meaningful. Performance degrades in strong mixture areas or when intra-class variability cannot be captured by a single reference spectrum scaled by per-pixel coefficients. As with other region-based methods, results depend on parameter settings and on the validity of the local noise approximation.

## Availability

CHOUCROUTE is released as open-source MATLAB code. The repository includes a demo script and documentation to reproduce a complete run from synthetic data.

## References
