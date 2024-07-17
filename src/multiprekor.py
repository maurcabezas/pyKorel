import argparse
import os
import numpy as np
import astropy.io.fits as fits
import astropy.constants as const
from astropy.io import ascii
from pathlib import Path
import astropy.units as u
from specutils.spectra import SpectralRegion
from specutils.manipulation import extract_region
from specutils import Spectrum1D
import configparser

class MultiPrekor:
    def __init__(self, prekorpath):
        self.prekorpath = f'{prekorpath}/etc/prekor.par'
        self.prekorpar = os.path.join(prekorpath, 'etc/prekor.par')
        # print 

    def load_input_data(self, input_file):
        return np.loadtxt(input_file, dtype=str)

    def write_config(self, speclist, wav_low, wav_up, nbin, hjd):
        config_content = f'''[Prekor config]
        speclist = {speclist}
        wav_low = {wav_low}
        wav_up = {wav_up}
        nbin = {nbin}
        hjd = {hjd}'''
        print (self.prekorpar)
        with open(self.prekorpar, 'w') as f:
            f.write(config_content)

    def run_prekor(self,prekorpath):
        os.system(f"python {prekorpath}/src/prekor.py")

    def process_spectra(self, input_file,prekorpath):
        datspec = self.load_input_data(input_file)

        for i, line in enumerate(datspec):
            speclist = str(datspec[i, 0])
            wav_low = int(datspec[i, 1])
            wav_up = int(datspec[i, 2])
            nbin = int(datspec[i, 3])
            hjd = str(datspec[i, 4])

            self.write_config(speclist, wav_low, wav_up, nbin, hjd)
            print (wav_low,wav_up,nbin,hjd)
            wav_ini = wav_low
            wav_fin = wav_up
            wav_cen = int((wav_ini + wav_fin) / 2)

            self.run_prekor(prekorpath)
            os.chdir(f"{wav_cen}_{nbin}/")
            # os.system("./run")
            os.chdir("../")

def main(args):
    prekorpath = '/home/cabezas/pyKorel'
    multi_prekor = MultiPrekor(prekorpath)
    multi_prekor.process_spectra(args.input_file,prekorpath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process input spectra for Korel analysis.")
    parser.add_argument("-l", "--input_file", help="Path to the input file containing spectra information.", required=True)
    args = parser.parse_args()
    main(args)
