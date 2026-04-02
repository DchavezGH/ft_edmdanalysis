## Extended dynamic mode decomposition (eDMD) over nonlinear dictionaries for neural data pipelines
David Chavez Huerta
## Introduction

This software provides a FieldTrip-integrated implementation of the Extended Dynamic Mode Decomposition (eDMD) algorithm for multi-channel neural data (MEG, EEG, and iEEG) analysis. eDMD is a data-driven method that generates an approximation to the Koopman operator over a set of observable functions. 

The frequency content of the neural data is deduced from the eDMD modes. Power associated with these frequencies is also calculated. By binning the frequency content into relevant neural bands (alpha, beta, gamma...), peak frequencies per band are calculated, as well as total power per band. Through the mode decomposition inherent to eDMD, the user can reconstruct the data dynamics from the initial snapshot. The function offers several possibilities for the set of observable functions: Random Fourier features (RFF), Hermite polynomials, polynomial basis, and the identity observable.

The package is meant to integrate with the user's FieldTrip pipeline, enabling them to incorporate nonlinear methods into existing workflows.

An overview of eDMD and its relation to the Koopman operator can be found in **[Matthew O. Williams' paper](https://link.springer.com/article/10.1007/s00332-015-9258-5)**. A manuscript on our own eDMD algorithm is still to be published. 

## Statement of need

Analyzing human EEG and MEG datasets often requires methods capable of capturing nonlinear, time-varying dynamics that traditional linear models cannot fully represent. Extended Dynamic Mode Decomposition (eDMD) provides a framework based on Koopman theory for learning interpretable models in a feature space, enabling the extraction of coherent temporal patterns from high-dimensional, noisy neural recordings. However, existing implementations generally require substantial customization, miss sufficient documentation, or lack integration with standard neuroscience toolboxes.

This software addresses this gap by providing a neuroscience-focused eDMD pipeline built around nonlinear dictionaries, enabling efficient approximation of nonlinear dynamics in high-dimensional EEG data. The eDMD implementation includes parameter choices tailored to the characteristics of EEG data, while still allowing users to adjust them easily for different analysis goals or needs
Concise documentation explains the effect of each parameter, supporting operation and customization. This makes the software immediately usable for cognitive-neuroscience data analysis.

The package aims to allow researchers to incorporate eDMD-based modeling directly into well established pre-processing and analysis workflows. It supports the substitution or extension of dictionaries, facilitating exploration of task-specific nonlinear observables or alternative basis functions. By providing an integrated, accessible, and extensible implementation, the software fills a methodological need for researchers seeking to integrate Koopman-based techniques into their own analyses.

## Installation

### Requirements

* MATLAB (R2021a or later recommended)

* FieldTrip toolbox (recent stable release recommended)

### Install FieldTrip

* Download FieldTrip from the official FieldTrip website: https://www.fieldtriptoolbox.org/

* Add it to your MATLAB path and initalize:

`ft_defaults`

* To verify the installation:

`ft_version`

### Install This Package

* Clone the repository using Git:

`git clone https://github.com/DchavezGH/ft_edmdanalysis`

* Or download the ZIP file and extract it

* Add the repository folder to your MATLAB path:

`addpath(genpath('path_to_repo'))`

* Then ft_edmdanalysis can be called as any other FieldTrip function.

### Reproducibility

* For deterministic behavior when using random Fourier features:

`cfg.seed = 1;`

* Tested Configuration

   * MATLAB R2023b

   * FieldTrip (2024 stable release)

* Users are encouraged to report compatibility issues via this repository 


NOTE:
The computational complexity of eDMD depends on dictionary size, data sampling, number of channels, and stacking depth. Users should adjust cfg.D, cfg.poly_degree, and cfg.nstacks carefully for large trials.

## Basic Usage

The function `ft_edmdanalysis` performs Extended Dynamic Mode Decomposition (eDMD) on FieldTrip raw data structures.  
The input data must be a valid FieldTrip **raw structure** containing time series data in `datain.trial`.

The function uses **`cfg.output`** to control the type of result produced.  
Three output modes are available:

- `'freq'` – interpolated spectrum derived from Koopman eigenvalues  (default)
- `'binned_peak_freq'` – power aggregated within frequency bins with peak frequency detection  
- `'raw'` – reconstructed time-domain signal using the Koopman model

---
### Using random Fourier features dictionary
```
cfg = [];
cfg.dictionary = 'rff';

dataout = ft_edmdanalysis(cfg, datain);
```
This runs eDMD using a random Fourier features dictionary on all trials in datain and returns a spectrum derived from Koopman eigenvalues. Random Fourier features are also the default dictionary.

### Using identity dictionary
```
cfg = [];
cfg.dictionary = 'identity';

dataout = ft_edmdanalysis(cfg, datain);
```
This runs classical DMD and returns a FieldTrip frequency structure containing the Koopman-derived spectrum.

### Using a Polynomial Dictionary
```
cfg = [];
cfg.dictionary = 'poly';
cfg.poly_degree = 3;

dataout = ft_edmdanalysis(cfg, datain);
```
This constructs a monomial polynomial dictionary up to degree 3.

### Using a Hermite Dictionary
```
cfg = [];
cfg.dictionary = 'hermite';
cfg.hermite_degree = 3;

dataout = ft_edmdanalysis(cfg, datain);
```
This applies probabilists’ Hermite polynomials independently to each channel.

## Function Output

As default, ft_edmdanalysis will output a FieldTrip freq structure, containing a spectrum generated from the Koopman decomposition. Details on the output structures are in the following section. Two other output options are available:  

### Binned Frequency Representation

To aggregate Koopman mode power into predefined frequency bands:
```
cfg = [];
cfg.output = 'binned_peak_freq';
cfg.freqEdges = [0 4 8 12 15 30 100];

dataout = ft_edmdanalysis(cfg, datain);
```

This returns a FieldTrip frequency structure where:

- power is summed within each frequency bin
- the dominant frequency per bin is stored separately

Peak frequencies per trial are stored in dataout.peakfreq

### State Reconstruction

To reconstruct the time-domain signal from Koopman modes:
```
cfg = [];
cfg.output = 'raw';

dataout = ft_edmdanalysis(cfg, datain);
```
This returns a FieldTrip raw structure containing reconstructed signals.

Optional normalization of the reconstruction can be enabled with `cfg.normalize_recon = true`

## Output Structures

The structure of `dataout` depends on `cfg.output`.


### 1. `cfg.output = 'freq'`

Returns a FieldTrip frequency structure representing the interpolated Koopman spectrum. This is the default output when `cfg.output` is empty (`[]`).

Fields:

* dataout.label = {'edmd'}
* dataout.freq = frequency vector (cfg.foi)
* dataout.powspctrm = [Ntrials x 1 x Nfreq]
* dataout.dimord = 'rpt_chan_freq'
* dataout.cfg

Additional diagnostics:

* dataout.modefreqs
* dataout.modepowers
* dataout.rank

Each entry corresponds to one trial.


### 2. `cfg.output = 'binned_peak_freq'`

Returns a FieldTrip frequency structure where power is aggregated in frequency bins.

Fields:

* dataout.label = {'edmd'}
* dataout.freq = bin centers
* dataout.powspctrm = [Ntrials x 1 x Nbins]
* dataout.dimord = 'rpt_chan_freq'
* dataout.cfg

Additional information:

* dataout.peakfreq


### 3. `cfg.output = 'raw'`

Returns a FieldTrip raw structure containing reconstructed signals.

Fields:

* dataout.label
* dataout.trial
* dataout.time
* dataout.fsample
* dataout.dimord
* dataout.cfg


## API

The main entry point of the package is:

`ft_edmdanalysis(cfg, datain)`


#### Inputs:

`cfg` – configuration structure controlling the decomposition

`datain` – FieldTrip raw data structure


#### Key Configuration Options


* Dictionary:  `cfg.dictionary`= `'identity'`, `'rff'`, `'poly'`, or `'hermite'`

*  Output selection: `cfg.output` = `'freq'`, `'binned_peak_freq'`, or `'raw'`
#### Dictionary parameters:

`cfg.D` number of random Fourier features
`cfg.gamma` RFF scaling parameter
`cfg.poly_degree` polynomial degree
`cfg.hermite_degree` Hermite polynomial degree

#### Algorithm parameters:

`cfg.nstacks` Hankel stacking depth  
`cfg.MA` moving average window 
`cfg.cut` rank truncation threshold
`cfg.seed` RNG seed for RFF reproducibility

#### Spectral parameters:

`cfg.foi` frequencies of interest
`cfg.freqEdges` frequency bin edges
`cfg.smooth` smoothing window for spectrum

#### Reconstruction normalization option:

`cfg.normalize_recon` true / false



### Community Guidelines

We welcome contributions, feedback, and questions from the community. Please be kind and constructive in all interactions

Contributing

* Contributions are welcome through pull requests

* Please fork the repository and create a new branch for your changes

* Ensure that your code is clear, documented, and consistent with the existing style

* Provide a clear description of the changes in your pull request

Reporting Issues

* Bugs and issues should be reported via the repository issue tracker

* Please include a clear description of the problem, relevant code or configuration
  
Support

* For questions or usage help, please open an issue in the repository. We will do our best to respond, but response times may vary


### License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).
See the LICENSE file for details.
