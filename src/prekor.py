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

class Prekor:
    def __init__(self, parpath):
        self.parpath = f'{parpath}/etc/prekor.par' #parpath
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
        config.read_file(open(f'{self.parpath}'))
        return config

    def snr_derived(self, spectrum, region=None):
        if region is not None:
            calc_spectrum = extract_region(spectrum, region)
        else:
            calc_spectrum = spectrum

        if hasattr(spectrum, 'mask') and spectrum.mask is not None:
            flux = calc_spectrum.flux[~calc_spectrum.mask]
        else:
            flux = calc_spectrum.flux

        n = len(flux)

        if n > 4:
            signal = np.median(flux)
            noise  = 0.6052697 * np.median(np.abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
            return signal / noise
        else:
            return 0.0

    def process_spectra(self):
        if Path("korel.dat").exists():
            os.rename("korel.dat", "old_korel.dat")

        print("Calculating SNR  -+-+-+-+")
        for i, spec_in_fits in enumerate(self.datspec):
            try:
                spec_in_ascii = self.convert_to_ascii(spec_in_fits)
                wav_in, flux_in = self.read_spectrum(spec_in_ascii)
                spec = Spectrum1D(spectral_axis=wav_in * u.AA, flux=flux_in * u.Unit('erg cm-2 s-1 AA-1'))
                limit = (self.wsnrl, self.wsnru)
                region = SpectralRegion(limit[0] * u.angstrom, limit[1] * u.angstrom)
                snr = self.snr_derived(spec, region=region)
                self.snrlist.append(snr)
            except Exception as e:
                print(f"Error processing {spec_in_fits}: {e}")

        if self.snrlist:
            snrmax = np.amax(self.snrlist)
        else:
            snrmax = 0

        print("SNR  -+-+-+-+ Done")
        print("\nSUMMARY:")
        fname = f"{self.wav_cen}_{self.nbin}"
        os.makedirs(f"{fname}/asc", exist_ok=True)
        os.makedirs(f"{fname}/model", exist_ok=True)

        prekor = []
        print("\nSpec name  HJD         SNR     Weight  Initial  final  wav.")
        rvstep = None  # Initialize rvstep to None

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
                print(f"{spec_in_ascii} {hjd:12.5f} {snr:7.3f} {weight:.3f}  {wav_out[0]:.3f} {wav_out[-1]:.3f}")
            except Exception as e:
                print(f"Error processing {spec_in_fits}: {e}")

        ascii.write([prekor], f"{fname}/tmp.res", overwrite=True, format='no_header')
        os.system(f"sed 's/\"//g' {fname}/tmp.res > {fname}/prekor.res")
        os.remove(f"{fname}/tmp.res")
        print("\nNumber of spectra: ", len(self.datspec))
        print("Max snr: ", snrmax)
        if rvstep is not None:
            print("RV step: ", rvstep)
        else:
            print("RV step: not calculated due to lack of valid spectra")
        if Path('korel.dat').exists():
            os.rename('korel.dat', f"{fname}/korel.dat")
        self.copy_files(fname)

    def convert_to_ascii(self, spec_in_fits):
        wav, flux, hdr, sp = self.read_spec(spec_in_fits)
        p = Path(spec_in_fits)
        spec_in_ascii = f"{p.stem}.asc"
        ascii.write([wav, flux], spec_in_ascii, overwrite=True, format='no_header')
        return spec_in_ascii

    def read_spec(self, filename):
        sp = fits.open(filename)
        hdr = sp[0].header
        flux = sp[0].data
        w = WCS(hdr, naxis=1, relax=False, fix=False)
        wave = w.wcs_pix2world(np.arange(len(flux)), 0)[0]
        return wave, flux, hdr, sp

    def get_hjd(self, spec_in_fits):
        with fits.open(spec_in_fits) as hdul:
            hjd = hdul[0].header.get(self.hjdhead)
        if hjd is None:
            raise ValueError(f"HJD header {self.hjdhead} not found in {spec_in_fits}")

        if (int(hjd/2400000) == 1):
            hjd = (float(hjd - 2400000))

        return hjd

    def read_spectrum(self, spec_in_ascii):
        try:
            wav, flux = np.loadtxt(spec_in_ascii, unpack=True)
        except Exception as e:
            raise ValueError(f"Error reading {spec_in_ascii}: {e}")
        return wav, flux

    def copy_files(self, fname):
        os.system(f"cp {self.parpath}/etc/korel.par {fname}/")
        os.system(f"cp {self.parpath}/etc/korel.gnu {fname}/")
        os.system(f"cp {self.parpath}/etc/rv.gnu {fname}/")
        os.system(f"cp {self.parpath}/etc/finalpaper.gnu {fname}/")

if __name__ == "__main__":
    parpath = '/home/cabezas/pyKorel'
    prekor = Prekor(parpath)
    prekor.process_spectra()
