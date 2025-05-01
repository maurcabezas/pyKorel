# pyKorel

`pyKorel` is a Python-based analysis tool designed to enhance the workflow of the original [Korel](https://www.astro.sk/~had/korel.html) software (Hadrava, 2004) for spectral disentangling of binary and multiple stellar systems. It provides a modern, user-friendly interface and additional functionalities without modifying the core Korel Fortran code. `pyKorel` can be considered a "Korel Suite," enabling a complete spectral analysis pipeline for disentangling composite spectra in Fourier space.

## Introduction to Korel

[Korel](https://www.astro.sk/~had/korel.html) is a powerful software package for disentangling the composite spectra of binary and multiple stellar systems. By performing spectral disentangling in Fourier space, Korel separates the contributions of individual components, allowing astronomers to study their properties (e.g., radial velocities, orbital parameters) in detail. However, the original Korel code relies on an MS-DOS-based interface, which can be cumbersome for modern workflows. `pyKorel` addresses this by providing Python tools to streamline data preparation, execution, and batch processing.

## Features

- **Prekor**: A Python module for preparing input spectra for Korel, improving efficiency over the original MS-DOS version. It handles FITS files, computes signal-to-noise ratios (SNR), and generates logarithmically resampled spectra.
- **Multiprekor**: A batch-processing tool for preparing multiple spectra files simultaneously, ideal for large datasets.
- **Korel.py**: A Python interface for running the Korel software, simplifying parameter setup and execution.
- **Multikorel**: A utility for managing multiple Korel runs, enabling parallel or sequential processing of different spectral regions or configurations.
- **Debug Mode**: Supports a debug mode (e.g., `python prekor.py -d`) for detailed logging, aiding in troubleshooting and development.

## Prerequisites

Before installing `pyKorel`, ensure you have the following:

- **Python**: Version 3.8 or higher (3.10 recommended).
- **Operating System**: Linux or macOS (Windows support is experimental due to Korel’s Fortran dependencies).
- **Dependencies**:
  - `numpy`
  - `matplotlib`
  - `astropy`
  - `specutils`
  - `configparser`
- **Fortran Compiler**: Required for compiling the original Korel code (e.g., `gfortran`).
- **Korel Source Code**: The original Korel Fortran code must be obtained separately (contact the author, [Petr Hadrava](https://www.astro.sk/~had/), for access).

## Installation

The easiest way to install `pyKorel` is by cloning the repository and running the provided `install.py` script. It is highly recommended to create an isolated Python virtual environment to avoid dependency conflicts.

### Step 1: Clone the Repository

```bash
git clone https://github.com/maurcabezas/pyKorel.git
cd pyKorel
```

### Step 2: Create a Virtual Environment

Create and activate a Python virtual environment to manage dependencies:
```bash
cd /path/to/your/preferred/location
python3.10 -m venv pykorel
source /path/to/your/preferred/location/pykorel/bin/activate
```
Note: You must activate the virtual environment every time you work with pyKorel:

```bash
source /path/to/your/preferred/location/pykorel/bin/activate
```

### Step 3: Install pyKorel

With the virtual environment activated, run the installation script:

```bash
python install.py
```

The install.py script will:

- Install Python dependencies (numpy, matplotlib, astropy, specutils, configparser).
- Install the pyKorel components: Prekor, Multiprekor, Korel.py, and Multikorel.
- Attempt to compile the Korel Fortran code (requires gfortran and the Korel source code in the appropriate directory).

Note: If you do not have the Korel source code, the installation will skip compiling Korel but still install the Python tools. You must manually place the compiled Korel executable in the appropriate directory (e.g., pyKorel/src/korel).

### Step 4: Testing Prekor

To verify that pyKorel is installed correctly, use the test data provided in the spec_test directory. The following steps test the Prekor module, which is the core component for preparing Korel input spectra.



#### Navigate to the spec_test directory:

```bash
cd spec_test
```

 #### Run the Prekor script:

```bash
python ../src/prekor.py
```
#### Expected Output:

- The script will process the test FITS files listed in speclist.txt (e.g., spec01_HJD2450000.00_combined.fits).
- It will compute the SNR for each spectrum, resample the spectra logarithmically, and generate output files in a directory named after the central wavelength and number of bins (e.g., 6562_1000).
- Console output will include:

```bash
INFO: Calculating SNR  -+-+-+-+
INFO: SNR  -+-+-+-+ Done
INFO: 
SUMMARY:
INFO: 
Spec name  HJD         SNR     Weight  Initial  final  wav.
INFO: spec01_HJD2450000.00_combined.asc 2450000.00000 150.000 0.750  6553.000 6570.000
INFO: spec10_HJD2450008.90_combined.asc 2450008.90000 199.500 1.000  6553.000 6570.000
INFO: Number of spectra: 10
INFO: Max snr: 199.5
INFO: RV step: 0.123
```


#### Output Files:
- 6562_1024/korel.dat: Korel input file with HJD, wavelengths, radial velocity steps, weights, and flux data.
- 6562_1024/prekor.res: Summary file with spectrum names, HJD, SNR, weights, and wavelength ranges.
- 6562_1024/asc/*.asc: Resampled ASCII spectra.
- 6562_1024/model/: Directory for model outputs (populated later by Korel).

Here 6562 is the central wavelenght and 1024 the number of pixels/bins.


### Debug Mode (Optional): 

To see detailed debugging output (e.g., wavelength ranges, flux values, SNR calculations), run:

```bash
python ../src/prekor.py -d
```

This will display additional information, such as:

```bash
DEBUG: FITS file: spec10_HJD2450008.90_combined.fits
DEBUG: WCS keywords: CRVAL1=6552.0, CDELT1=0.010005002501250625, CRPIX1=1, CUNIT1=Angstrom
DEBUG: Raw WCS wavelength range: 6.552e-07 to 6.572e-07 (WCS output)
DEBUG: Warning: WCS wavelengths are too small (6.572e-07). Converting from meters to Angstrom.
DEBUG: Corrected wavelength range: 6552.0 to 6572.0 Angstrom
DEBUG: Writing ASCII file: spec10_HJD2450008.90_combined.asc
DEBUG: Wavelength range before writing: 6552.0 to 6572.0 Angstrom
DEBUG: Reading ASCII file: spec10_HJD2450008.90_combined.asc
DEBUG: Wavelength range from ASCII: 6552.0 to 6572.0 Angstrom
DEBUG: Flux range: 0.310937801492317 to 1.013250519305709
DEBUG: Spectrum spec10_HJD2450008.90_combined.fits wavelength range: 6552.0 to 6572.0 Angstrom
DEBUG: Input spectrum: Spectrum1D (length=2000)
DEBUG: Region: Spectral Region, 1 sub-regions: (6555.0 Angstrom, 6568.0 Angstrom)
DEBUG: Extracted calc_spectrum: Spectrum1D (length=1300)
DEBUG: Extracted flux: [0.99901126 0.9906245 ...] erg / (Angstrom s cm2)
DEBUG: Flux shape: (1300,), Flux values: [0.99901126 0.9906245 ...]
DEBUG: Signal: 0.9862870173012894, Noise: 0.004942383381187716
```

Troubleshooting:
- If you see errors (e.g., missing FITS files, invalid HJD headers), check that the speclist.txt file lists valid FITS files and that the prekor.par configuration file is correctly set up.
- If the SNR is zero for any spectrum, run with -d to inspect the debug output for issues like empty flux arrays or wavelength mismatches.
- Ensure the FITS files have valid WCS headers (e.g., CRVAL1, CDELT1, CUNIT1) and HJD keywords.

### Usage
#### Configuration

Edit the etc/prekor.par file to specify parameters for Prekor:

```bash
[Prekor config]
speclist = speclist.txt
wav_low = 6550
wav_up = 6575
nbin = 1000
hjd = HJD
```

- speclist: Path to a text file listing FITS files (one per line).
- wav_low, wav_up: Wavelength range in Angstroms (e.g., 6550–6575 Å for the H-alpha region).
- nbin: Number of bins for logarithmic resampling (e.g., 1000).
- hjd: FITS header keyword for Heliocentric Julian Date (e.g., HJD).

For Korel.py, configure the korel.par file in the output directory (e.g., 6562_1000/korel.par) with parameters like orbital elements and disentangling options. Refer to the Korel documentation for details.

## FAQ

- Q: What if I don’t have the Korel source code? A: pyKorel’s Python tools (Prekor, Multiprekor) will still work, but you need the compiled Korel executable to run Korel.py or Multikorel. Contact Petr Hadrava to obtain Korel.

- Q: Why is the SNR zero for some spectra? A: This could indicate a wavelength mismatch, empty flux arrays, or invalid FITS headers. Run with -d to inspect debug output and check prekor.par and FITS headers.

- Q: How do I configure korel.par? A: Refer to the Korel documentation and the sample korel.par in etc/. Adjust parameters like orbital period, eccentricity, and wavelength range based on your system.


## Citation

If you use pyKorel in your research, please cite:

- **Hadrava, P. (2004)**. "Fourier Disentangling of Composite Spectra." Publications of the Astronomical Institute of the Czech Academy of Sciences, 92, 15.
- **Cabezas, M. (2025)**. pyKorel: A Python Suite for Korel Spectral Disentangling. GitHub: https://github.com/maurcabezas/pyKorel.
