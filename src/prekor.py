import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from pathlib import Path
import astropy.units as u
from specutils.spectra import SpectralRegion, Spectrum1D
from specutils.manipulation import extract_region
import configparser
from astropy.wcs import WCS
import logging
import argparse

class Prekor:
    def __init__(self, initpath):
        self.initpath = f'{initpath}/etc/prekor.par'
        self.config = self.load_config()
        self.speclist = self.config.get('Prekor config', 'speclist')
        self.wav_low = float(self.config.get('Prekor config', 'wav_low'))
        self.wav_up = float(self.config.get('Prekor config', 'wav_up'))
        self.nbin = int(self.config.get('Prekor config', 'nbin'))
        self.hjdhead = str(self.config.get('Prekor config', 'hjd'))
        self.wav_ini = self.wav_low
        self.wav_fin = self.wav_up
        self.wav_cen = int((self.wav_ini + self.wav_fin) / 2)
        self.dwav = self.wav_fin - self.wav_ini
        self.wsnrl = self.wav_low + 2
        self.wsnru = self.wav_up - 2
        self.snrlist = []
        self.datspec = np.loadtxt(self.speclist, dtype='str')
        self.nspec = len(self.datspec)

    def load_config(self):
        config = configparser.ConfigParser()
        config.read_file(open(f'{self.initpath}'))
        return config

    def snr_derived(self, spectrum, region=None):
        logging.debug(f"Input spectrum: {spectrum}")
        logging.debug(f"Region: {region}")
        if region is not None:
            calc_spectrum = extract_region(spectrum, region)
            logging.debug(f"Extracted calc_spectrum: {calc_spectrum}")
            logging.debug(f"Extracted flux: {calc_spectrum.flux}")
            logging.debug(f"Extracted spectral axis: {calc_spectrum.spectral_axis}")
        else:
            calc_spectrum = spectrum
            logging.debug(f"Using full spectrum: {calc_spectrum}")

        flux = calc_spectrum.flux
        if hasattr(calc_spectrum, 'mask') and calc_spectrum.mask is not None:
            flux = flux[~calc_spectrum.mask]
            logging.debug(f"Applied mask, flux shape: {flux.shape}")

        logging.debug(f"Flux shape: {flux.shape}, Flux values: {flux}")
        if len(flux) == 0:
            logging.error("Flux array is empty")
            return 0.0
        if np.all(np.isnan(flux)) or np.all(flux == 0):
            logging.error("Flux contains only NaNs or zeros")
            return 0.0

        n = len(flux)
        if n > 4:
            signal = np.nanmedian(flux)
            noise_array = np.abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n])
            noise = 0.6052697 * np.nanmedian(noise_array)
            logging.debug(f"Signal: {signal}, Noise: {noise}")
            if noise == 0 or np.isnan(noise) or np.isnan(signal):
                logging.error("Invalid signal or noise value")
                return 0.0
            return signal / noise
        else:
            logging.error(f"Flux length {n} is too short")
            return 0.0

    def process_spectra(self):
        if Path("korel.dat").exists():
            os.rename("korel.dat", "old_korel.dat")

        logging.info("Calculating SNR  -+-+-+-+")
        for i, spec_in_fits in enumerate(self.datspec):
            try:
                spec_in_ascii = self.convert_to_ascii(spec_in_fits)
                hjd = self.get_hjd(spec_in_fits)
                wav_in, flux_in = self.read_spectrum(spec_in_ascii)
                spec = Spectrum1D(spectral_axis=wav_in * u.AA, flux=flux_in * u.Unit('erg cm-2 s-1 AA-1'))
                logging.debug(f"Spectrum {spec_in_fits} wavelength range: {wav_in.min()} to {wav_in.max()} Angstrom")

                # Validate SpectralRegion
                limit = (self.wsnrl, self.wsnru)
                if not (wav_in.min() <= self.wsnru and wav_in.max() >= self.wsnrl):
                    logging.error(f"SpectralRegion {limit} does not overlap with spectrum range {wav_in.min()} to {wav_in.max()}")
                    self.snrlist.append(0.0)
                    continue

                region = SpectralRegion(limit[0] * u.angstrom, limit[1] * u.angstrom)
                snr = self.snr_derived(spec, region=region)
                self.snrlist.append(snr)
            except Exception as e:
                logging.error(f"Error processing {spec_in_fits}: {e}")

        if self.snrlist:
            snrmax = np.amax(self.snrlist)
        else:
            snrmax = 0

        logging.info("SNR  -+-+-+-+ Done")
        logging.info("\nSUMMARY:")
        fname = f"{self.wav_cen}_{self.nbin}"
        os.makedirs(f"{fname}/asc", exist_ok=True)
        os.makedirs(f"{fname}/model", exist_ok=True)

        prekor = []
        logging.info("\nSpec name  HJD         SNR     Weight  Initial  final  wav.")
        rvstep = None

        for i, spec_in_fits in enumerate(self.datspec):
            try:
                spec_in_ascii = self.convert_to_ascii(spec_in_fits)
                hjd = self.get_hjd(spec_in_fits)
                wav_in, flux_in = self.read_spectrum(spec_in_ascii)
                spec = Spectrum1D(spectral_axis=wav_in * u.AA, flux=flux_in * u.Unit('erg cm-2 s-1 AA-1'))
                region = SpectralRegion(limit[0] * u.angstrom, limit[1] * u.angstrom)
                snr = self.snr_derived(spec, region=region)
                weight = snr / snrmax if snrmax > 0 else 0

                wav_new = np.linspace(self.wav_ini, self.wav_fin, num=self.nbin)
                wav_log = 3e5 * np.log(wav_new / wav_new[0])
                dwavlog = [wav_log[j + 1] - wav_log[j] for j in range(self.nbin - 1)]
                d_wav_log = np.average(dwavlog)
                rvstep = d_wav_log
                equidistant_log_scale = np.array([wav_log[0] + j * d_wav_log for j in range(self.nbin)])
                wav_out = np.exp(equidistant_log_scale / 3e5) * wav_new[0]
                flux_out = np.interp(wav_out, wav_in, flux_in)

                os.remove(spec_in_ascii)
                ascii.write([wav_out, flux_out], f"{fname}/asc/{spec_in_ascii}", overwrite=True, format='no_header')

                with open('korel.dat', 'a') as koreldat:
                    if i != 0:
                        koreldat.write('\n')
                    koreldat.write(f'{hjd:12.5f}{wav_out[0]:10.4f}{rvstep:7.3f}  {weight:5.3f}     {self.nbin}\n')
                    for j, flux in enumerate(flux_out):
                        if j % 10 == 0 and j != 0:
                            koreldat.write('\n')
                        koreldat.write(f' {flux:7.5f}')

                prekor.append(f'{spec_in_ascii} {hjd:12.5f} {snr:7.3f} {weight:.3f} {wav_out[0]:.3f} {wav_out[-1]:.3f}')
                logging.info(f"{spec_in_ascii} {hjd:12.5f} {snr:7.3f} {weight:.3f}  {wav_out[0]:.3f} {wav_out[-1]:.3f}")
            except Exception as e:
                logging.error(f"Error processing {spec_in_fits}: {e}")

        ascii.write([prekor], f"{fname}/tmp.res", overwrite=True, format='no_header')
        os.system(f"sed 's/\"//g' {fname}/tmp.res > {fname}/prekor.res")
        os.remove(f"{fname}/tmp.res")
        logging.info(f"\nNumber of spectra: {len(self.datspec)}")
        logging.info(f"Max snr: {snrmax}")
        if rvstep is not None:
            logging.info(f"RV step: {rvstep}")
        else:
            logging.info("RV step: not calculated due to lack of valid spectra")
        if Path('korel.dat').exists():
            os.rename('korel.dat', f"{fname}/korel.dat")
        self.copy_files(fname)

    def convert_to_ascii(self, spec_in_fits):
        wav, flux, hdr, sp = self.read_spec(spec_in_fits)
        p = Path(spec_in_fits)
        spec_in_ascii = f"{p.stem}.asc"
        logging.debug(f"Writing ASCII file: {spec_in_ascii}")
        logging.debug(f"Wavelength range before writing: {wav.min()} to {wav.max()} Angstrom")
        ascii.write([wav, flux], spec_in_ascii, overwrite=True, format='no_header')
        return spec_in_ascii

    def read_spec(self, filename):
        sp = fits.open(filename)
        hdr = sp[0].header
        flux = sp[0].data
        w = WCS(hdr, naxis=1, relax=False, fix=False)
        pixel_indices = np.arange(len(flux))
        wave = w.wcs_pix2world(pixel_indices, 0)[0]
        
        logging.debug(f"FITS file: {filename}")
        logging.debug(f"WCS keywords: CRVAL1={hdr.get('CRVAL1')}, CDELT1={hdr.get('CDELT1')}, CRPIX1={hdr.get('CRPIX1')}, CUNIT1={hdr.get('CUNIT1')}")
        logging.debug(f"Raw WCS wavelength range: {wave.min()} to {wave.max()} (WCS output)")
        
        if hdr.get('CUNIT1') == 'Angstrom' and wave.max() < 1e-3:
            logging.debug(f"Warning: WCS wavelengths are too small ({wave.max()}). Converting from meters to Angstrom.")
            wave *= 1e10
        
        logging.debug(f"Corrected wavelength range: {wave.min()} to {wave.max()} Angstrom")
        return wave, flux, hdr, sp

    def get_hjd(self, spec_in_fits):
        with fits.open(spec_in_fits) as hdul:
            hjd = hdul[0].header.get(self.hjdhead)
        if hjd is None:
            raise ValueError(f"HJD header {self.hjdhead} not found in {spec_in_fits}")

        if int(hjd / 2400000) == 1:
            hjd = float(hjd - 2400000)
        return hjd

    def read_spectrum(self, spec_in_ascii):
        try:
            wav, flux = np.loadtxt(spec_in_ascii, unpack=True)
            logging.debug(f"Reading ASCII file: {spec_in_ascii}")
            logging.debug(f"Wavelength range from ASCII: {wav.min()} to {wav.max()} Angstrom")
            logging.debug(f"Flux range: {flux.min()} to {flux.max()}")
            return wav, flux
        except Exception as e:
            raise ValueError(f"Error reading {spec_in_ascii}: {e}")

    def copy_files(self, fname):
        os.system(f"cp {self.initpath}/etc/korel.par {fname}/")
        os.system(f"cp {self.initpath}/etc/korel.gnu {fname}/")
        os.system(f"cp {self.initpath}/etc/rv.gnu {fname}/")
        os.system(f"cp {self.initpath}/etc/finalpaper.gnu {fname}/")

def setup_logging(debug=False):
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format='%(levelname)s: %(message)s'
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prekor: Process spectra for Korel disentangling')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug output')
    args = parser.parse_args()
    setup_logging(args.debug)
    
    initpath = 'update/path/installation'
    prekor = Prekor(initpath)
    prekor.process_spectra()
