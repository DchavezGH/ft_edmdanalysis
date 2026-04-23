---
title: 'Extended dynamic mode decomposition over nonlinear dictionaries integrated in FieldTrip for neurophysiological data pipelines'
tags:
  - FieldTrip
  - MATLAB
  - Koopman operator
  - eDMD
  - DMD
  - Random Fourier features
  - EEG
  - MEG
authors:
  - name: David Chavez-Huerta
    orcid: 0000-0003-2054-7751
    affiliation: 1
    corresponding: true 
  - name: Mohammad Khosravi
    orcid: 0000-0002-4873-1115
    affiliation: 1
    corresponding: false
affiliations:
 - name: Delft University of Technology, Netherlands
   index: 1
   ror: 02e2c7k09 
   date: 23 april 2026
bibliography: paper.bib
license: GPL-3.0
---


# Summary
This software provides a MATLAB implementation of the Extended Dynamic Mode Decomposition (eDMD) algorithm integrated with the FieldTrip software toolbox [@Oostenveld2010FieldTripOS] for neural data analysis. The presented software offers several sets of observable functions: Random Fourier features (RFFs), Hermite sequences, a polynomial basis, and the identity observable.
Through the mode decomposition inherent to eDMD, the user can generate a reconstruction of the dynamics in the data. 
Frequency content is deduced from the eDMD modes. Power associated to these frequencies is also calculated. By binning frequency content into standard brainwave bands (i.e., delta, theta, alpha, beta, gamma), peak frequencies per band are calculated, as well as total power per band. The specific outputs (a signal reconstruction, a spectrum, or a binned dominant-frequency) are generated as standard Fieldtrip data structures so eDMD can be immediately integrated in existing analysis pipelines.   

The package integrates with FieldTrip, enabling users to incorporate nonlinear methods into existing EEG workflows that currently rely on linear algorithms. The software also allows for the implementation of custom dictionaries, making it adaptable to task-specific and model-specific nonlinear representations. The software is meant to support researchers studying neural dynamics, cognitive processes, and data-driven modeling in neuroscience. 


# Statement of need

Analyzing human EEG and MEG datasets often requires methods capable of capturing nonlinear, time-varying dynamics that traditional linear models cannot fully represent. Extended Dynamic Mode Decomposition (eDMD) provides a framework based on Koopman theory for learning interpretable models in a feature space [@Williams2014ADA], enabling the extraction of coherent temporal patterns from high-dimensional, noisy neural recordings. However, existing implementations of the eDMD algorithm generally require substantial customization, lack adequate documentation, or are not integrated with standard neuroscience toolboxes.

This software addresses this gap by providing a neuroscience-focused eDMD pipeline built around nonlinear dictionaries, enabling efficient approximation of nonlinear dynamics in high-dimensional EEG data. The eDMD implementation includes parameter choices tailored to the characteristics of EEG data, while still allowing users to adjust them easily for different analysis goals or applications. Concise documentation explains the effect of each parameter, supporting operation and customization. This makes the software readily usable for cognitive-neuroscience data analysis.

The package is integrated with FieldTrip, allowing researchers to incorporate eDMD-based modeling directly into well-established pre-processing and analysis workflows. It supports the substitution or extension of dictionaries, facilitating exploration of task-specific nonlinear observables or alternative basis functions. By providing a compatible, accessible, and extensible implementation, the software fills a methodological need for researchers seeking to integrate Koopman-based techniques into their own analyses.  


# State of the Field

Decomposition methods based on the Koopman operator (such as eDMD) have become popular analysis tools in several fields [@budivsic2012applied]. In the last two decades, the inherent potential of the operator (to describe the behavior of complex dynamical systems in terms of a linear operator) has motivated many researchers to propose algorithms and to explore applications [@brunton2021modern]. In particular, high-dimensional systems with nonlinear dynamics in the fields of climate science, fluid dynamics, neuroscience, and neural-network training have been investigated with promising results [@ghosh2024koopman]. This has resulted in plenty of software related to Koopman-based methodologies.

PyKoopman [@Pan2024] is an extensive Python package for computing data-driven approximations to the Koopman operator. It implements multiple dictionaries and several methods (including DMD, DMDc, eDMD, and kDMD). The well-known KoopmanLab [@xiong2023koopmanlab] is a Python implementation focused on solving partial differential equations (PDE) by approximating the Koopman operator with a Krylov sequence. A deep learning approach to learning Koopman eigenfunctions [@lusch2018deep] is available in the DeepKoopman repository, also written in Python. There are also some related MATLAB implementations, such as an eDMD toolbox based on Orthogonal Polynomials [@GarciaTenorio2022AMT], and an implementation of the measure-preserving version of eDMD described in [@Colbrook2022TheMA]. In addition to these established options, there exists a considerable number of independent projects that address specific problems with Koopman-based methods. 

For our own research purposes, there was a need for a MATLAB-based alternative easily integrable into neural data pipelines (ideally, integrated with the FieldTrip toolbox). We decided to build an eDMD script, capable of switching between dictionary options, designed with trial-based electrophysiology data in mind. In this context, dynamics reconstruction is important, but the focus is on the extraction of frequency content (associated with the Koopman modes) aligned with standard EEG analysis. The main ambition of the code is to promote the use of nonlinear dictionaries in eDMD as an alternative to standard methods (PCA, ICA, FFT) in neuroscience analysis pipelines by providing an eDMD tool integrated in FieldTrip. To the best of our knowledge, there is no such alternative available at the moment.


# Software Design

Neurophysiological recordings (MEG, EEG, fMRI, ECoG) are typically organized as trial-based, multichannel time series, which introduces specific requirements for data handling, reproducibility, noise rejection, and interpretability. Our implementation therefore operates directly on FieldTrip raw data structures and their associated configuration structures, processing each trial independently while preserving experiment structure, data provenance, and compatibility with existing pre-processing pipelines. In addition to computing Koopman decompositions, the software produces frequency features and band-power derived from the Koopman eigenvalues. These are generated as standard FieldTrip data structures, aligning the output with common analysis practices in neuroscience. 

A central design choice is the use of a configurable dictionary, that is, a set of observables, for constructing the lifted state space required by eDMD. The implementation supports multiple dictionary types, including identity (equivalent to classical DMD), RFFs, polynomial bases, and Hermite polynomials. This allows users to tailor the balance between expressiveness, interpretability, and computational cost. For example, RFFs provide scalable nonlinear embeddings for high-dimensional data, while polynomial and Hermite dictionaries offer more structured and interpretable representations at the expense of increased dimensionality. By exposing the dictionary choice and its parameters through a unified configuration interface, the software enables systematic exploration of Koopman models across different regimes without modifying the underlying algorithm. Default options are carefully chosen to support users who are less familiar with Koopman analysis: a RFF dictionary with a reasonable cardinality for trials in the range of 10-20 seconds, and a conservative stacking parameter that supports good transient reproduction without excessive computational expense.

There is more than one way to implement eDMD: The algorithm in our code is designed to ensure stability and computational efficiency in the high-dimensional settings typical of neural data. A compact singular value decomposition (SVD) formulation is first used to reduce the data rank, then the Koopman operator can be computed without explicitly constructing large matrices. In the current implementation, 70\% of the total running time is dedicated to SVD calculation, while the remaining operations have been streamlined to limit additional computational overhead.
 
Additional design choices, such as optional data smoothing, reproducible random feature generation, and normalization within specific dictionaries, further enhance flexibility. The software also prioritizes practical usability by directly using the Koopman eigenvalues to build a frequency-domain representation. These design choices aim to balance theoretical exploration of the Koopman operator with the practical constraints of analyzing large-scale neural series.

# Research impact statement
The software is currently used to identify theta and high-beta activity in the frontal midline cortex area correlated with working memory encoding in an EEG dataset via the RFF dictionary. These features are subsequently used in downstream statistical analyses to quantify the correlation between spectral activity and behavioral measures. Specifically, the extracted Koopman-based features are used to assess the presence of high-beta band activity during memory encoding events. Preliminary analyses indicate an increase in total power in the high-beta band associated to memory capacity (that is, the more items in memory, the more high-beta power) [@chavez_edmd_rff_working_memory]. This demonstrates the utility of the software as a feature-extraction software tool within standard neuroscience analysis pipelines.

To facilitate adoption by the community, the toolbox is distributed as an independent open-source repository and it is documented on the official website of the FieldTrip toolbox [@Oostenveld], where it is presented in the Extensions section. The integration into the FieldTrip ecosystem increases the visibility and accessibility of the method to an already existing user base and demonstrates interoperability with established neuroscience analysis pipelines.


# AI usage disclosure

Microsoft Copilot was used on occasion during the development of this work to troubleshoot isolated coding issues, brainstorm alternative approaches, and explore potential submission venues prior to submission. All AI‑generated suggestions were critically reviewed, tested for technical accuracy, and incorporated only when they aligned with the authors’ own criteria. The authors have access to Microsoft Copilot enterprise data protection, preventing their input from being used to train the model. Nonetheless, the authors remain fully responsible for the submitted manuscript and software.

# Acknowledgements

We acknowledge the contributions of Bart De Schutter in the review of this manuscript. Likewise, we acknowledge the assistance with integration and compatibility of Jan-Mathijs Schoffelen and Robert Oostenveld from the FieldTrip developing team. We also acknowledge the support of the Delft Center for Systems and Control, Faculty of Mechanical Engineering, TU Delft.

# References