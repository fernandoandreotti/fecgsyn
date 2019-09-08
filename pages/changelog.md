---
layout: default
title: Changelog
priority: 4
---

## Changelog

Here is a list of changes made on each version of the <em>FECGSYN</em> simulator.

- **v1.3-pre-release : 2019.09.05**<br>
Asymmetric volume conductor modeling code added by @emersonkeenan
  - **Features**
    - Import 3D anatomic models as presented in Keenan et al 2018
    - Set vernix caseosa distribution and thickness in predefined regions
    - Generate finite element mesh from 3D models using Iso2Mesh with compartment specific mesh size and conductivity
    - Calculate lead field matrices using FieldTrip-Simbio pipeline
    - Simulate ECG signals at chosen sensor positions

- **v1.2 : 2017.03.13** <br>
  Update to include novel fetal signal quality metrics
  - **Features**
    - Novel signal quality indices (SQI) added as presented in Andreotti _et al_ 2017
    - FHR correction using SQI and Kalman filtering added as presented in Andreotti _et al_ 2017
    - Naive Bayes classifier summarizes different SQI metrics into a consensus metric.
    - 4 additional QRS detectors modified for FECG use
    - Minor restructuring of folders

- **v1.1 : 2016.03.15** <br>
  Major update and inclusion of several extraction methods and  benchmarking tools.
  - **Features**
    - Benchmark toolbox for FQRS detection and morphological analysis
    - 8 Extraction methods included (as described in Andreotti et al 2016)
    - FECGSYNDB released
    - Overall restructuration of toolbox directory
  - **Bug Fixes**
    - 0.1 : @fernandoandreotti fixed several bugs in calibration, generation and evaluation    
    - 0.3 : noise generation function updated by Julien Oster
    - 0.4 : @mosalvi GUI updated to match new release?

- **v1.0 : 2014.09.15** <br>
  First major versioned release by @jbehar and @fernandoandreotti
  - **Features**
    - Each dipole modelled as an individual source
    - Simulator included translation and rotation of cardiac dipoles
    - Realistic AR-modelled noise sources added
    - Calibration routine for differen SNR levels
    - Example containing patho-physiological events
    - @mohsalvi developed GUI for easy use    
