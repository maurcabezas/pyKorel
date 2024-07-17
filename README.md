# pyKorel

pyKorel is an analysis tool designed to work with the original Korel code (Hadrava, 2004). It aims to provide a more user-friendly interface and additional functionalities without altering the core Korel source code.

## Introduction to Korel

Korel is a powerful software package for disentangling composite spectra of binary and multiple stellar systems. It performs spectral disentangling in the Fourier space, which is crucial for studying such systems' components individually.

## Features

- **Prekor**: A complete Python method for preparing Korel input spectra, improving the workflow over the original MS-DOS based version.
- **Multiprekor**: A batch processing tool for handling multiple spectra files.
- **Korel.py**: The main interface for running the Korel software.
- **Multikorel**: A utility for managing multiple Korel runs.

## Installation

Clone the repository and install the dependencies:

```sh
git clone https://github.com/yourusername/py-Korel.git
cd py-Korel
pip install -r requirements.txt
