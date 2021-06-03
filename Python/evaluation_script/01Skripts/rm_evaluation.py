# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Maximilian Sender
# Copyright   : Copyright 2021, Institute of Chemical Engineering, Prof. Dr. Dirk Ziegenbalg, Ulm University'
# License     : GNU LGPL
# =============================================================================
"""The module (rm = radiometric measurement) rm_evaluation contains the class RM_Evaluation, which helps evaluation
radiometric measurements performed with a radiometric scanning device developed at the Institute of Chemical
Engineering, University Ulm. By creating an instance of the class, it can lead you through the procedure.

Run your python instance in the Folder "01Scripts" and put th measurement files th the folder "02Data".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
# =============================================================================
# Imports
# =============================================================================
import os
import pprint
import logging
import time
import pandas as pd
from pathlib import Path
import re
import math
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler
import numpy as np
import scipy
from scipy.optimize import minimize
from supplementary.colour_system import cs_hdtv as cs
from lmfit.lineshapes import skewed_gaussian
import sympy as sp

logging.basicConfig(format="%(asctime)s: %(message)s", level=logging.INFO, datefmt="%H:%M:%S")


class RM_Evaluation:
    """
    The class RM_Evaluation can be used for further evaluation of measurements performed by the scanning device.
    By initializing an instance of the class the script will scan existing subfolders for exported RM-data. This Data
    must be stored in a specific Folder structure: (Normally created by the export function of the measurement script.)

    folder: 02Data
    |
    |-- folder: date + time + reactor name
    |   |
    |   |-- folder: z-position 1
    |   |   |-- file: *.feather ....................... contains measured spectra in mW nm**-1
    |   |   |-- file: *.log ........................... contains all set parameters for measurement
    |   |   |-- heatmap pictures and *.csv file ....... not used by this script
    |   |
    |   |-- folder: z-position 2
    |   |   |-- ...
    |   |
    |   |-- ...
    |
    |-- ...

        -------------------
        Attributes constructed from parameters by __init___:

            saveram (boolean, optional, default = False):
                If True, photon flux spectra are not loaded to use less RAM.

            wl_ntsr_borders (tuple of two integers, optional):
                For the calculation of the noise to signal ratio the spectrum gets separated in two regimes. By default
                the regime used for defining the noise are the spectrometer pixels 0-150. The signal regime are the pixels 150 to -1
                (end of list). This default is defined by the method import_raw. Other borders should be defined here as
                tuple of integers if the automatic import routine is used. The tuple defines the pixel separating the
                two regimes (def: 150) and the end of the signal regime.

            predef (list of strings, optional):
                Defines a predefined answer to the import questions of the class when initialized. You give a list of
                eight the answers, it will be red from right to left (list.pop).

                When initializing an instance of the class and predef == False you will be asked
                if you want to choose a Directory to evaluate:

                Choose a Directory to evaluate?
                >?|
                
                Two answers are possible:
                • When answering with ’y’ , all export folders in the same directory which have the correct structure and
                content will be listed. The program then will ask for the number of the folder
                which should be imported.
                • If ’n’ is given as answer, all importable .feather files in all subdirectories are listed. The
                files can be picked individually by their number.
                After choosing the desired data, the program will ask whether it should import these files:
                
                import raw?
                >?|
                
                If the answer is ’n’ , you can start the import procedure manually with custom parameters by calling
                the function import_raw . Otherwise the program will also ask if colours should be processed:
                
                process colors?
                >?|
                
                and for a value of the ntsr_border :
                keep ntsr border of 0.550? Or type new
                >?|

                After that, all files will be red and processed. Note that this can afford some RAM. The answers to the
                import questions can be also given in advance by the predef parameter with the class initialization.
                In case you want to import the second folder ( answer2 = 2 ) in the directory ( answer1 = ’y’ )
                and import( answer3 = ’y’ ) without processing the colours ( answer4 = ’n’ ), while keeping the
                ntsr_border unchanged ( answer5 = ’y’ ), the list [’y’, ’n’, ’y’, 2, ’y’])must be passed
                as parameter predef . The program reads this list from right to left.

        -------------------
        Output Attributes:

            alldfs (numpy array,
                    shape: (number of z positions, number of x-pixels, number of y-pixels, number of spectrometer pixels),
                    dtype='float32'):
                Array containing the measured spectra in mW nm**-1.

            all_photon_dfs (numpy array,
                    shape: (number of z positions, number of x-pixels, number of y-pixels, number of spectrometer pixels),
                    dtype='float32'):
                Array containing the measured spectra in photons nm**-1 s**-1.


            allwaves (numpy array,
                      shape: (number of z positions, number of spectrometer pixels),
                      dtype='float32'):
                Array containing the respective wavelengths in nm for the spectrometer pixels.

            all_vars (numpy array,
                      shape: (number of z positions),
                      dtype='float32'):
                Array containing a dict with the respective measurement variables:
                Name': Name of the Scan
                'z_pos': z-position
                'Size': Size of the scan im pixels.
                int_time': Integration time.

            all_integrals (numpy array,
                           shape: (number of z positions, number of x-pixels, number of y-pixels),
                           dtype='float32'):
                Containing the received power per pixel in mW.

            all_photon_integrals (numpy array,
                                  shape: (number of z positions, number of x-pixels, number of y-pixels),
                                  dtype='float32'):
                Containing the received number of photons per pixel.

            all_rec_power (numpy array,
                          shape: (number of z positions),
                          dtype='float32'):
                Containing the received power per scan in mW.

            all_rec_photons (numpy array,
                             shape: (number of z positions),
                             dtype='float32'):
                Containing the received number of photons per scan.

            all_nev_abs (numpy array,
                         shape: (number of z positions, number of x-pixels, number of y-pixels, number of spectrometer
                         pixels),
                         dtype='float32'):
                Array containing the measured spectra of not absorbable power in mW nm**-1.

            all_nev_abs_photons (numpy array,
                                 shape: (number of z positions, number of x-pixels, number of y-pixels, number of
                                 spectrometer pixels),
                                 dtype='float32'):
                Array containing the measured spectra of not absorbable photons.

            all_ratios (numpy array,
                        shape: (number of z positions, number of x-pixels, number of y-pixels),
                        dtype='float32'):
                Array containing the peak ratio for every measurement pixel.

            all_colors (numpy array,
                        shape: (number of z positions, number of x-pixels, number of y-pixels, red values, green values,
                        blue values),
                        dtype='float32'):
                Array containing the RGB values.

            all_integrals_nev_abs (numpy array,
                                   shape: (number of z positions, number of x-pixels, number of y-pixels),
                                   dtype='float32'):
                Containing the not absorbable received power per pixel in mW.

            all_integrals_remaining (numpy array,
                                     shape: (number of z positions, number of x-pixels, number of y-pixels),
                                     dtype='float32'):
                Containing the not absorbed received power per pixel in mW.

            all_photon_integrals_remaining (numpy array,
                                            shape: (number of z positions, number of x-pixels, number of y-pixels),
                                            dtype='float32'):
                Containing the not absorbed received photons per pixel.

    """
    def __init__(self, saveram=False, wl_ntsr_borders=False, predef=False):
        self.script_starting_time = time.time()
        self.saveram = saveram
        self.dir_path = os.path.dirname(os.path.realpath(__file__))+'\\..\\02Data'
        self.foundthefiles = []
        self.thefiles = []
        self.dirs = []
        self.thedir = ''
        self.frames_red = False
        self.ntsr_border = 0.55
        self.waveleghth_borders_fo_ntsr = wl_ntsr_borders
        self.cal_int_time = 0.03551
        self.planck_constant = 6.626 * 10 ** (-34)
        self.speed_of_light = 2.998 * 10 ** 8
        self.avogadro_number = 6.022 * 10 ** 23
        self.alldfs = np.array([], dtype='float32')
        self.all_photon_dfs = np.array([], dtype='float32')
        self.allwaves = np.array([], dtype='float32')
        self.all_vars = np.array([])
        self.all_integrals = np.array([], dtype='float32')
        self.all_photon_integrals = np.array([], dtype='float32')
        self.all_rec_power = np.array([], dtype='float32')
        self.all_rec_photons = np.array([], dtype='float32')
        self.all_nev_abs = np.array([], dtype='float32')
        self.all_nev_abs_photons = np.array([], dtype='float32')
        self.all_ratios = np.array([], dtype='float32')
        self.all_colors = np.array([])
        self.df_for_int_exp = pd.DataFrame({})
        self.color_cycle = cycler('color', ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple',
                                            'pink', 'brown', 'orange', 'teal', 'coral', 'lightblue', 'lime', 'lavender',
                                            'turquoise', 'darkgreen', 'tan', 'salmon', 'gold', 'orchid', 'crimson',
                                            'darkblue'])
        self.color_cycle_list = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple', 'pink',
                                 'brown', 'orange', 'teal', 'coral', 'lightblue', 'lime', 'lavender', 'turquoise',
                                 'darkgreen', 'tan', 'salmon', 'gold', 'orchid', 'crimson', 'darkblue']
        #  Here we scan all folders for feather data and let the user choose which ones to import
        for self.root, dirs, files in os.walk(self.dir_path):
            for file in files:
                if file.endswith('.feather'):
                    self.foundthefiles.append(str(self.root+'\\'+str(file)))
                    if Path(self.root).parent in self.dirs:
                        pass
                    else:
                        self.dirs.append(Path(self.root).parent)

        if len(self.foundthefiles) == 1:
            pass
            self.thefiles = self.foundthefiles
        else:
            print('Choose a Directory to evaluate?')
            if predef is not False:
                bool_dir = predef.pop()
            else:
                bool_dir = input()
            logging.info('%s' % bool_dir)
            if bool_dir == 'y':
                self.chose_dir = True
                print('I found these Folders:\n')
                counter = 1
                for i in self.dirs:
                    print(str(counter) + ':   ' + i.parts[-1])
                    counter += 1
                print('\nWhich one should be processed? Only choose one.')
                if predef is not False:
                    self.user = predef.pop()
                else:
                    self.user = int(input())
                print('you chose the path: %s' % self.dirs[self.user-1])
                self.thedir = self.dirs[self.user-1]
                for i in self.foundthefiles:
                    if Path(i).parts[0:-2] == Path(self.dirs[self.user-1]).parts[:]:
                        print('Adding %s' % Path(i).parts[-1])
                        self.thefiles.append(i)
            else:
                self.chose_dir = False
                print('I found these Files:\n')
                counter = 1
                for i in self.foundthefiles:
                    print(str(counter) + ':   ' + i)
                    counter += 1
                print(str(counter) + ':   all')
                print('\nWhich ones should be processed? Separate with blanks.')
                if predef is not False:
                    self.user = predef.pop()
                else:
                    self.user = input()
                    self.user = self.user.split()
                for i in range(0, len(self.user)):
                    self.user[i] = int(self.user[i])
                if counter in self.user:
                    self.thefiles = self.foundthefiles
                else:
                    for i in self.user:
                        self.thefiles.append(self.foundthefiles[i-1])
        if predef is not False:
            import_raw = predef.pop()
        else:
            import_raw = input('import raw?')
        if import_raw == 'y':
            if predef is not False:
                process_colors = predef.pop()
            else:
                process_colors = input('process colors?')
            if process_colors == 'y':
                input_colors = True
            else:
                input_colors = False
            if predef is not False:
                input_ntsr = predef.pop()
            else:
                input_ntsr = input('keep ntsr border of %.3f? Or type new' % self.ntsr_border)
            try:
                self.ntsr_border = float(input_ntsr)
            except ValueError:
                logging.info('ntsr border not changed')
            self.import_raw(colors=input_colors)

    def normalize(self, arr):
        """Helper method to normalize an array."""
        arr_min = np.min(arr)
        return (arr-arr_min)/(np.max(arr)-arr_min)

    def find_nearest(self, array, value):
        """Helper method, finds the index of a e.g. wavelength"""
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def pos_to_num(self, pos):
        """Helper method, calculates a position on a heatmap in the respective pixel."""
        return int(round(pos/0.39))

    def cos_corr_point(self, array, center, distance):
        """
        Returns an array containing cosine correction values for a
        point light source, e.g. LED
        ------------------
        Parameters:

            array (2D numpy array):
                The array which should be corrected, here only the correct shape is necessary. Using the actual array
                is useful considering the ’auto’ option of the parameter center.

            center (tuple or 'auto':
                The center of the light source on the array. If center is ’auto’ , the position of the largest value in
                the array array is used. Otherwise, pass a tuple of coordinates in cm.

            distance (number):
                The distance between the measurement canvas and the light source in cm.
        """
        if 'auto' in center:
            center = np.reshape(np.array(np.where(array == np.amax(array)))*0.39+(0.39/2), 2)
        self.olsh_positions = np.indices(array.shape)*0.39+(0.39/2)
        self.positions = np.zeros((array.shape[0], array.shape[1], 2))
        self.positions[:, :, 0] = self.olsh_positions[0, :, :]
        self.positions[:, :, 1] = self.olsh_positions[1, :, :]  # we need to give the positions array a different shape
        self.vectors = np.subtract(self.positions, center)  # we calculate the vectors pointing from the lightsource to the pixel
        print('vectors shape: ' + str(self.vectors.shape))
        self.norms = np.linalg.norm(self.vectors, axis=-1)
        print('norms shape: ' + str(self.norms.shape))
        return 1/np.cos(np.arctan(self.norms/distance))

    def cos_corr_stick(self, array, x_center=70, y_min=84, y_max=107, distance=4.5):
        """
        Returns an array containing cosine correction values for a longish light source. See cos_corr_point.
        ------------------
        Parameters:

            array (2D numpy array):
                The array which should be corrected, here only the correct shape is necessary. Using the actual array
                is useful considering the ’auto’ option of the parameter center.

            x_center (integer, default = 70):
                Pixel position of a longish light source in x direction.

            y_min (integer, default = 84):
                Lower pixel position of a longish light source in y direction.

            y_max (integer, default = 107):
                Upper pixel position of a longish light source in y direction.

            distance (number):
                The distance between the measurement canvas and the light source in cm.
        """
        x_center = x_center*0.39+(0.39/2)
        y_min = y_min*0.39+(0.39/2)
        y_max = y_max*0.39+(0.39/2)
        olsh_positions = np.indices(array.shape)*0.39+(0.39/2)
        positions = np.zeros((array.shape[0], array.shape[1], 2))
        positions[:, :, 0] = olsh_positions[0, :, :]
        positions[:, :, 1] = olsh_positions[1, :, :]

        y_max_helper = np.zeros(positions.shape)
        y_max_helper[..., 0] = y_max_helper[..., 1] = np.where(positions[:, :, 0] >= y_max, 1, 0)

        y_min_helper = np.zeros(positions.shape)
        y_min_helper[..., 0] = y_min_helper[..., 1] = np.where(positions[:, :, 0] <= y_min, 1, 0)

        y_mid_helper = np.zeros(positions.shape)
        y_mid_helper[..., 0] = y_mid_helper[..., 1] = np.where((positions[:, :, 0] > y_min) & (positions[:, :, 0] < y_max), 1, 0)
        print(y_mid_helper)
        vectors_y_max = np.subtract(positions, [y_max, x_center])*y_max_helper
        vectors_y_min = np.subtract(positions, [y_min, x_center])*y_min_helper
        mid_positions = np.ones(positions.shape)
        mid_positions[:, :, 0] = mid_positions[:, :, 0]*y_max
        mid_positions[:, :, 1] = positions[:, :, 1]
        vectors_y_mid = np.subtract(mid_positions, [y_max, x_center])*y_mid_helper
        vectors = np.where(vectors_y_max == 0, vectors_y_min, vectors_y_max)
        vectors = np.where(vectors == 0, vectors_y_mid, vectors)
        print('vectors shape: ' + str(vectors.shape))
        norms = np.linalg.norm(vectors, axis=-1)
        print('norms shape: ' + str(norms.shape))
        return 1/np.cos(np.arctan(norms/distance))

    def vars_from_log(self, file):
        """
        Helper method to extract variables from the logfile.
        """
        import re
        filepath = file
        search_int_time = 'int_time'
        matched_int_time = ''

        search_Name = 'Name:'
        matched_name = ''

        search_Size = 'Size:'
        matched_size = ''

        search_Z = 'z_pos:'
        matched_z = ''

        with open(filepath, 'r') as file:
            for line in file:
                if search_int_time in line:
                    matched_int_time = float(re.sub("\D+[_]\D+", "", line.rstrip()))

        with open(filepath, 'r') as file:
            for line in file:
                if search_Name in line:
                    matched_name = re.sub("\D+[: ]", "", line.rstrip())

        with open(filepath, 'r') as file:
            for line in file:
                if search_Size in line:
                    matched_size = int(re.sub("\D", "", line.rstrip()))

        with open(filepath, 'r') as file:
            for line in file:
                if search_Z in line:
                    matched_z = float(re.sub("\D+[_]\D+", "", line.rstrip()))
        out = {'int_time': matched_int_time, 'Name': matched_name, 'Size': matched_size, 'z_pos': matched_z}
        pprint.pprint(out)
        return out

    def noisetosignal(self, a, borders, axis=-1):
        """
        Helper method for the import_raw method, calculates the noise to signal ratio for an array.
        ------------
        Parameters:
            a (array, shape: (no. of. scans, no. of x-pixels, no. of y-pixels, spec- pixels):
                The array.
            borders (tuple, shape:(2)):
                The borders for ntds calculation, see docstring of this class.
            axis (integer, default = -1):
            Axis parameter for the numpy std method, no need to change.
        """
        c = np.std(a[:, :, :, :borders[0]], axis=axis)
        d = np.std(a[:, :, :, borders[0]:borders[1]], axis=axis)
        return np.where(d == 0, 1, c/d)

    def show_histogram(self, values):
        """
        Shows a histogram of the given array.
        """
        n, bins, patches = plt.hist(values.reshape(-1), 50, density=1)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])

        for c, p in zip(self.normalize(bin_centers), patches):
            plt.setp(p, 'facecolor', cm.viridis(c))

        plt.show()

    def normalize(self, arr):
        """
        Helper method for show_histogram.
        """
        arr_min = np.min(arr)
        return (arr-arr_min)/(np.max(arr)-arr_min)

    def SI_plot(self, comment=''):
        """
        Plots and saves a figure of supporting information, subdivided in 4 parts. Subfigures 1,2 and 3 display
        histograms of the NTSR, peak ratios and P of all scans. The latter subfigure displays the received power per
        scan over the z-position and the integration time over the z-position. The received power over the z-position
        is fitted with a linear regression and the function is given.
        """
        fig, axes = plt.subplots(ncols=2,
                                 nrows=2,
                                 figsize=(8, 6),
                                 sharex=False,
                                 sharey=False)
        """overall Title"""
        fig.suptitle(self.all_vars[0]['Name']+' SI plot '+comment)

        """Making the Histogram of all ntsr"""
        target = (0, 0)
        axes[target].set_prop_cycle(self.color_cycle)
        axes[target].tick_params(labelsize=7)
        axes[target].set_title('histogram of signal to noise ratios', fontsize=9)
        axes[target].set_xlabel("signal to noise / 1", fontsize=9, labelpad=0.1)
        axes[target].set_ylabel("frequency", fontsize=9, labelpad=0.1)
        if self.separate_borders_mode is True:
            axes_twin_ntsr = axes[target].twinx()
            axes_twin_ntsr.set_prop_cycle(self.color_cycle)
        for i in range(0, len(self.thefiles)):
            n, bins, patches = axes[target].hist(self.all_ntsr[i].reshape(-1), 50, density=1, histtype='step',
                                                 label=str(self.all_vars[i]['z_pos'])+'$\\,$cm')
            if self.separate_borders_mode is True:
                axes_twin_ntsr.plot([self.separate_ntsr[i, 0, 0], self.separate_ntsr[i, 0, 0]], [0, n.max()/2], '-.',
                                    label='%.3f' % self.separate_ntsr[i, 0, 0], alpha=0.5)
                axes_twin_ntsr.set_yticks([])
        if self.separate_borders_mode is False:
            axes[target].plot([self.ntsr_border, self.ntsr_border], [0, n.max()], 'k-.'  # , label='border'
                              )
            xmin, xmax = axes[target].get_xlim()
            ymin, ymax = axes[target].get_ylim()
            t_xpos = xmax+0.01*(xmax-xmin)
            t_ypos = ymin+0.01*(ymax-ymin)
            axes[target].text(t_xpos, t_ypos, 'ntsr_border = %.3f' % self.ntsr_border, rotation='vertical', fontsize=6)
        if self.separate_borders_mode is True:
            axes_twin_ntsr.legend(fontsize=4.5, ncol=2, labelspacing=0.01, loc='center left', title='ntsr borders', title_fontsize=4.5, columnspacing=1)

        """Making the Histogram of all ratios"""
        target = (0, 1)
        axes[target].tick_params(labelsize=7)
        axes[target].set_prop_cycle(self.color_cycle)
        axes[target].set_title('histogram of peak ratios', fontsize=9)
        axes[target].set_xlabel("peak ratio/ 1", fontsize=9, labelpad=0.1)
        axes[target].set_ylabel("frequency", fontsize=9, labelpad=0.1)
        for i in range(0, len(self.thefiles)):
            n, bins, patches = axes[target].hist(self.all_ratios[i].reshape(-1), 50, density=1, histtype='step',
                                                 label=str(self.all_vars[i]['z_pos'])+'$\\,$cm', log=True)
        axes[target].legend(fontsize=4.5, ncol=2, labelspacing=0.01, loc='upper right', columnspacing=1)

        """plotting power over distance"""
        target = (1, 0)
        axes[target].tick_params(labelsize=7)
        axes[target].set_prop_cycle(self.color_cycle)
        axes[target].set_title('received power and\nintegration time over z', fontsize=9)
        axes[target].set_xlabel("z / cm", fontsize=9, labelpad=0.1)
        axes[target].set_ylabel("P / W", fontsize=9, labelpad=0.1)
        axes_twin = axes[target].twinx()
        axes_twin.tick_params(labelsize=7)
        axes_twin.set_ylabel("integration time / ms", fontsize=9, labelpad=0.1)
        x = np.zeros(len(self.thefiles))
        int_times = np.zeros(len(self.thefiles))
        # marker: https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/marker_reference.html
        marker_style = dict(linestyle=':', color='0.8', markersize=7, mfc="C0", mec="C0")
        marker_style2 = dict(linestyle=':', color='0.8', markersize=7, mfc="C1", mec="C1")
        for i in range(0, len(self.thefiles)):
            x[i] = self.all_vars[i]['z_pos']
            int_times[i] = self.all_vars[i]['int_time']*1000
        rec = axes[target].plot(x, self.all_rec_power*1E-6, label='received power', marker='.', **marker_style)
        abs = axes[target].plot(x, self.all_nev_abs*1E-6, label='never absorbed power', marker=7, **marker_style)
        int = axes_twin.plot(x, int_times, label='integration time', marker='2', **marker_style2)
        # legend: https://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
        all = rec+abs+int
        labs = [l.get_label() for l in all]

        """fitting the slope of data"""
        coef_z = np.polyfit(x, self.all_rec_power*1E-6, 1)
        coef_z_na = np.polyfit(x, self.all_nev_abs*1E-6, 1)
        coef_int = np.polyfit(int_times, self.all_rec_power*1E-6, 1)
        coef_int_na = np.polyfit(int_times, self.all_nev_abs*1E-6, 1)
        poly1d_fn_z = np.poly1d(coef_z)
        poly1d_fn_z_na = np.poly1d(coef_z_na)
        poly1d_fn_int = np.poly1d(coef_int)
        poly1d_fn_int_na = np.poly1d(coef_int_na)
        # poly1d_fn is now a function which takes in x and returns an estimate for y

        power_at_cal_time = poly1d_fn_int(self.cal_int_time)
        power_at_cal_time_na = poly1d_fn_int_na(self.cal_int_time)

        axes[target].plot(x, poly1d_fn_z(x), '--k')
        axes[target].plot(x, poly1d_fn_z_na(x), '--k')
        xmin, xmax = axes[target].get_xlim()
        ymin, ymax = axes[target].get_ylim()
        t_xpos = xmin+0.5*(xmax-xmin)  # the relative position on the xaxis
        t_ypos = ymin+0.8*(ymax-ymin)  # the relative position on the yaxis
        axes[target].text(t_xpos, t_ypos, '$m_{rec}$ = %.2f$\\cdot 10^{-3}\\,$Wcm$^{-1}$\nor %.2f$\\cdot 10^{-3}\\,$Wms$^{-1}$'
                                          '\n$P_{@caltime} = %.3f$' % (poly1d_fn_z[1]*1E3, poly1d_fn_int[1]*1E3,
                                                                       power_at_cal_time), fontsize=4.5)
        t_xpos = xmin+0.5*(xmax-xmin)  # the relative position on the xaxis
        t_ypos = ymin+0.3*(ymax-ymin)  # the relative position on the yaxis
        axes[target].text(t_xpos, t_ypos, '$m_{nev abs}$ = %.2f$\\cdot 10^{-3}\\,$Wcm$^{-1}$\nor %.2f$\\cdot 10^{-3}\\,$Wms$^{-1}$'
                                          '\n$P_{@caltime} = %.3f$' % (poly1d_fn_z_na[1]*1E3, poly1d_fn_int_na[1]*1E3,
                                                                       power_at_cal_time_na), fontsize=4.5)

        axes[target].legend(all, labs, fontsize=4.5, ncol=1, labelspacing=1, loc='best',
                            bbox_to_anchor=(0.05, 0.25, 0.5, 0.5), columnspacing=1)

        """Making the Histogram of all integrals"""
        target = (1, 1)
        axes[target].tick_params(labelsize=7)
        axes[target].set_prop_cycle(self.color_cycle)
        axes[target].set_title('histogram of receives power\nper pixel', fontsize=9)
        axes[target].set_xlabel("P / \u03BCW ", fontsize=9, labelpad=0.1)
        axes[target].set_ylabel("frequency", fontsize=9, labelpad=0.1)
        for i in range(0, len(self.thefiles)):
            n, bins, patches = axes[target].hist(self.all_integrals[i].reshape(-1), 50, density=1, histtype='step',
                                                 label=str(self.all_vars[i]['z_pos'])+'$\\,$cm', log=True)
        axes[target].legend(fontsize=4.5, ncol=2, labelspacing=0.01, loc='upper right', columnspacing=1)

        fig.subplots_adjust(
                # left=0.0,  # the left side of the subplots of the figure
                # right=0.9,  # the right side of the subplots of the figure
                bottom=0.1,  # the bottom of the subplots of the figure
                top=0.9,     # the top of the subplots of the figure
                wspace=0.35,  # the amount of width reserved for space between subplots,
                              # expressed as a fraction of the average axis width
                hspace=0.4)

        if self.thedir == '':
            self.thedir = Path(self.dir_path)
        fig.savefig(self.thedir/Path(self.all_vars[0]['Name'] + 'SI_plot' + comment + '.png'),
                    dpi=300,
                    bbox_inches='tight'
                    )

    def read_frames(self):
        """
        Helper method for import_raw.
        """
        for counter, thefile in np.ndenumerate(self.thefiles):
            self.vars = self.vars_from_log(thefile.replace('feather', 'log'))
            if counter[0] == 0:
                self.all_vars = np.array([self.vars]*(len(self.thefiles)))
            self.all_vars[counter] = self.vars
        self.all_z = np.zeros(len(self.all_vars))
        for num, i in np.ndenumerate(self.all_vars):
            self.all_z[num] = self.all_vars[num]['z_pos']
        self.all_vars = np.take_along_axis(self.all_vars, np.argsort(self.all_z), axis=0)
        self.thefiles = np.take_along_axis(np.array(self.thefiles), np.argsort(self.all_z), axis=0)

        for counter, thefile in np.ndenumerate(self.thefiles):
            print(thefile)

            logging.info('started read')
            imported_df = pd.read_feather(thefile)
            logging.info('ended read')

            if counter[0] == 0:
                self.all_integrals = np.zeros((len(self.thefiles), self.vars['Size'], self.vars['Size']), dtype='float32')  # now already a 3d array, so the reshaped array sould be written to it
                self.all_ratios = np.zeros((len(self.thefiles), self.vars['Size'], self.vars['Size']), dtype='float32')  # now already a 3d array, so the reshaped array sould be written to it
                self.all_ntsr = np.zeros((len(self.thefiles), self.vars['Size'], self.vars['Size']), dtype='float32')  # now already a 3d array, so the reshaped array sould be written to it
                self.all_ntsr_helper = np.zeros((len(self.thefiles), self.vars['Size'], self.vars['Size']), dtype='float32')  # now already a 3d array, so the reshaped array sould be written to it
                self.all_colors = np.zeros((len(self.thefiles), self.vars['Size'], self.vars['Size'], 3), dtype='float32')  # now already a 3d array, so the reshaped array sould be written to it
                self.alldfs = np.zeros((len(self.thefiles), self.vars['Size'], self.vars['Size'], len(imported_df['wavelength'])), dtype='float32')  # now already a 3d array, so the reshaped array sould be written to it
                self.all_photon_dfs = np.zeros((len(self.thefiles), self.vars['Size'], self.vars['Size'], len(imported_df['wavelength'])), dtype='float32')  # now already a 3d array, so the reshaped array sould be written to it
                self.allwaves = np.zeros((len(self.thefiles), len(imported_df['wavelength'])), dtype='float32')

            try:
                new_trans = np.loadtxt(Path(".\supplementary\TransferCurves\Transfercurve.csv"), delimiter=';')*(self.cal_int_time / self.all_vars[counter]['int_time'])
            except:
                raise Exception("Can\'t find the transfer curve")

            if np.array_equal(new_trans, imported_df['transfer_curve']) is True:
                new_data_df = imported_df.iloc[:, 3:]
            else:
                new_data_df = imported_df.iloc[:, 3:].multiply(new_trans/imported_df['transfer_curve'], axis=0)

            self.alldfs[counter] = np.reshape(np.transpose(new_data_df.values), (self.all_vars[counter]['Size'], self.all_vars[counter]['Size'], len(imported_df['wavelength'])))
            self.alldfs[counter][1::2] = np.flip(self.alldfs[counter][1::2], 1)
            self.allwaves[counter] = imported_df['wavelength']
            del imported_df
            del new_data_df
        self.alldfs = np.where(np.isnan(self.alldfs), 0.0, self.alldfs)
        self.frames_red = True

    def import_raw(self,
                   ntsr_border=None,
                   colors=False,
                   waveleghth_borders_fo_ntsr=(150, -1),
                   ntsr_helper_border=0,
                   separate_borders_mode=False,
                   analysis=False):
        """
        Starts the import procedure, normally automatically started within the class initialization.
        ------------------
        Parameters:

            ntsr_border (number, optional):
                Pixels with a NTSR (stds) smaller than this border will be set to 0 or black. If not set the value from
                class initialization is used.

            colors (boolean, optional, default is False):
                If True colours will be processed.

            waveleghth_borders_fo_ntsr (tuple of 2 integers, default is (150, -1):
                See parameter docstring for wl_ntsr_borders of this class.

            separate_borders_mode (boolean, optional, default is False):
                If True for every slice a new wavelength NTSR border is found. Namely the biggest stds which belongs to
                an integral value smaller than stds_helper_border.

            ntsr_helper_border (number, optional, default is 0):
                Needed for separate_borders_mode.

            analysis (boolean, optinal, default is False):
                If True, a log file is written showing the total received power over the z-positions.
        """
        if ntsr_border is None:
            ntsr_border = self.ntsr_border
        self.ntsr_border = ntsr_border
        self.separate_borders_mode = separate_borders_mode
        self.ntsr_helper_border = ntsr_helper_border
        if self.frames_red is False:
            self.read_frames()
        if self.waveleghth_borders_fo_ntsr is False:
            self.waveleghth_borders_fo_ntsr = waveleghth_borders_fo_ntsr

        logging.info('start ntsr')
        self.all_ntsr = self.noisetosignal(self.alldfs, borders=self.waveleghth_borders_fo_ntsr, axis=-1)
        logging.info('end ntsr')

        logging.info('start photons')
        """making the photons in nanomol per second and nm"""
        print('ntsr_border = %.3f'%ntsr_border)
        self.conv_to_photons_array = self.allwaves[0]*1e-9*1e6/self.planck_constant/self.speed_of_light/self.avogadro_number
        self.all_photon_dfs = self.alldfs*self.conv_to_photons_array
        logging.info('end photons')

        logging.info('start integrating')
        self.all_integrals = np.trapz(self.alldfs, self.allwaves[0], axis=-1) * 0.1521
        logging.info('power done')
        self.all_photon_integrals = np.trapz(self.all_photon_dfs, self.allwaves[0], axis=-1) * 0.1521
        logging.info('photons done')
        if self.saveram is True:
            del self.all_photon_dfs
        self.all_ntsr_helper = np.where(self.all_integrals < ntsr_helper_border, self.all_ntsr, 0)

        """separate borders mode means, that for every slice a new ntsr border is found. Namely the biggest ntsr which
        belongs to a integral value smaller than ntsr_helper_border"""
        if separate_borders_mode is True:
            """sparate_ntsr is an array with the length of the numper of slices. Containing the biggest ntsr
            which belongs to an integral value smaller than ntsr_helper_border"""
            self.ntsr_border = 'separate'
            self.separate_ntsr = np.max(self.all_ntsr_helper, (1, 2))
            logging.info('separate ntsr borders mode, found this:' + str(self.separate_ntsr))
            self.separate_ntsr = np.ones(self.all_ntsr.shape)*np.reshape(self.separate_ntsr, (len(self.thefiles), 1, 1))
            self.all_integrals = np.where(self.all_ntsr < self.separate_ntsr, self.all_integrals, 0)
            self.all_photon_integrals = np.where(self.all_ntsr < self.separate_ntsr, self.all_photon_integrals, 0)

        if separate_borders_mode is False:
            self.all_integrals = np.where(self.all_ntsr < ntsr_border, self.all_integrals, 0)
            self.all_photon_integrals = np.where(self.all_ntsr < ntsr_border, self.all_photon_integrals, 0)

        if self.saveram is False:
            self.alldfs_nev_abs = self.alldfs/(10**(7*skewed_gaussian(self.allwaves[0], 16386626.6, 567.358141, 37193531.6, -2830490.97)))
            self.all_photon_dfs_nev_abs = self.alldfs_nev_abs*self.conv_to_photons_array
            self.all_photon_integrals_nev_abs = np.trapz(self.all_photon_dfs_nev_abs, self.allwaves[0], axis=-1) * 0.1521
            if separate_borders_mode is False:
                self.all_photon_integrals_nev_abs = np.where(self.all_ntsr < ntsr_border, self.all_photon_integrals_nev_abs, 0)
            if separate_borders_mode is True:
                self.all_photon_integrals_nev_abs = np.where(self.all_ntsr < self.separate_ntsr, self.all_photon_integrals_nev_abs, 0)
            self.all_integrals_nev_abs = np.trapz(self.alldfs_nev_abs, self.allwaves[0], axis=-1) * 0.1521
            if separate_borders_mode is False:
                self.all_integrals_nev_abs = np.where(self.all_ntsr < ntsr_border, self.all_integrals_nev_abs, 0)
            if separate_borders_mode is True:
                self.all_integrals_nev_abs = np.where(self.all_ntsr < self.separate_ntsr, self.all_integrals_nev_abs, 0)
            self.all_integrals_remaining = self.all_integrals-self.all_integrals_nev_abs
            self.all_photon_integrals_remaining = self.all_photon_integrals-self.all_photon_integrals_nev_abs
        logging.info('end integrating')
        logging.info('start making ratios')
        if separate_borders_mode is True:
            self.all_ratios = self.make_peak_difference(self.alldfs, ntsr_border=self.separate_ntsr, ntsr=self.all_ntsr)
        if separate_borders_mode is False:
            self.all_ratios = self.make_peak_difference(self.alldfs, ntsr_border=ntsr_border, ntsr=self.all_ntsr)
        logging.info('end making ratios')

        if colors is True:
            logging.info('start making colors')
            lam = np.arange(380., 781., 5)  # lambda table for spec_to_xyz
            for index, i in np.ndenumerate(np.zeros((len(self.thefiles), self.all_vars[0]['Size'], self.all_vars[0]['Size']))):
                self.all_colors[index] = cs.spec_to_rgb(np.interp(lam, self.allwaves[0], self.alldfs[index]))
            if separate_borders_mode is True:
                self.all_colors = self.all_colors * np.reshape(np.where(self.all_ntsr<self.separate_ntsr, 1, 0),(len(self.thefiles), self.all_vars[0]['Size'], self.all_vars[0]['Size'],1))
            if separate_borders_mode is False:
                self.all_colors = self.all_colors * np.reshape(np.where(self.all_ntsr<ntsr_border, 1, 0),(len(self.thefiles), self.all_vars[0]['Size'], self.all_vars[0]['Size'],1))
            logging.info('end making colors')

        if self.chose_dir is True:
            logging.info('chose_dir is True')
            filename = str(self.thedir/Path(time.strftime("%y%m%d_%H%M%S")+'_'+self.vars['Name'] +'_analysis.log'))
            outtxt = 'Pos / cm; Received Power / W; Not Absorbable Power / W\n'
            self.all_rec_power = np.zeros(len(self.thefiles))
            self.all_rec_photons = np.zeros(len(self.thefiles))
            self.all_nev_abs = np.zeros(len(self.thefiles))
            rec_sum = 0
            nev_abs_sum = 0
            for i in range(0, len(self.thefiles)):
                logging.info('writing %.3f to %i'%(np.sum(self.all_integrals[i]),i))
                self.all_rec_power[i] = np.sum(self.all_integrals[i])
                self.all_rec_photons[i] = np.sum(self.all_photon_integrals[i])
                if self.saveram is False:
                    self.all_nev_abs[i] = np.sum(self.all_integrals_nev_abs[i])
                outtxt += (str(self.all_vars[i]['z_pos'])+ '; ' + str(np.round(self.all_rec_power[i]*1E-6, 3)) + '; ' + str(np.round(self.all_nev_abs[i]*1E-6, 3))+'\n')
            outtxt += ('mean; '+str(np.round(np.sum(self.all_rec_power)/len(self.thefiles)*1E-6, 3))+ '; ' +str(np.round(np.sum(self.all_nev_abs)/len(self.thefiles)*1E-6, 3)))
            if analysis is True:
                with open(filename, 'w') as f:
                    f.write(outtxt)

    def mak_peak_difference_from_cache(self, dataframe, peak1, peak2, ntsr_border=False, total_width=2):
        """
        Can be used if you want to calculate a new peak difference after import.
        -------------------
        Parameters:

            dataframe (array of scans): The array containing the scans, eg. alldfs.

            peak1 (integer): Index of the first peak.

            peak2 (integer): Index of the second peak.

            ntsr_border (number, default is ntsr_border of class): The border to determine no signal.

            total_width( number, default = 2): Half width of peak heigt determination.
        """
        if ntsr_border is False:
            ntsr_border = self.ntsr_border
        self.all_ratios = self.make_peak_difference(dataframe, self.all_ntsr, ntsr_border, total_width, [peak1, peak2])

    def make_peak_difference(self, dataframe, ntsr, ntsr_border, width=2, peakindices=[274, 526]):
        """Helper method fpr import_raw, return an array of the ratio of the two peaks, 0 when no signal"""
        ratio = (np.mean(dataframe[:, :, :, peakindices[0]-width:peakindices[0]+width], axis=-1)/np.mean(dataframe[:, :, :, peakindices[1]-width:peakindices[1]+width], axis=-1))*np.where(ntsr < ntsr_border, 1, 0)
        ratio = np.clip(ratio, -10, 10)
        real_clips = np.array([-2, 8])
        ratio = np.clip(ratio, real_clips[0], real_clips[1])
        return ratio

    def plot(self, mode='integrals', data=False, comment = '', cbar_title='', bar_to_max=False):
        """
        Plots and saves heatmaps of all loaded scans.
        -------------------
        Parameters:

            mode (string, default is 'integrals'):
                Defines the plotted data and display mode (matplotlib), combinations of data and method are
                automatically recognized. Possible data is:

                ’integrals’
                    Heatmap of the received power per pixel.
                    Can be combined with ’nev_abs’ and ’remaining’.

                ’ratios’
                    Heatmap of the peak ratios per pixel.

                ’custom’
                    Heatmap of a custom 2D array, given by the parameter data.

                ’photons’
                    Heatmap of the received photons per pixel.
                    Can be combined with ’nev_abs’ and ’remaining’.

                ’color’
                    Image in heatmap style, showing the calculated color of the pixels.
                    !!! Can’t be combined with other output options!!!!

                ’integrals’ and ’photons’Can be combined with:

                    ’nev_abs’
                        Not absorbable photons or power.
                    ’remaining’.
                        Not absorbed photons or power.

                Except ’color’ all other options also can be combined with:

                    ’contour_clabel’
                        Contour heatmap with labels. See Matplotlib documetatation.

                    ’contourf’
                        Contourf heatmap. See Matplotlib documetatation.

                    ’contour’
                        Contour heatmap. See Matplotlib documetatation.

            data (array, default is False):
                Numpy array for mode == 'custom'.

            comment (string, default = ''):
                Additional string to be displayed in the figure.

            cbar_title (string, default = ''):
                Title of the colorbar for mode == 'custom'.

            bar_to_max (boolean, default is False):
                If True, all colorbars are scaled to the absolute maximum of the figure set.
        """
        boundsstring = ''
        """Here we set the size of the figure, when smaller than 6 we make two cols. when bogger, we make 3"""
        if len(self.thefiles) <= 6:
            cols = 2
        else:
            cols = 3
        fig, axes = plt.subplots(ncols=cols,
                                 nrows=math.ceil(len(self.thefiles)/cols),
                                 figsize=(cols*3.4, math.ceil(len(self.thefiles)/cols)*3),
                                 sharex=False,
                                 sharey=False)
        if cols < len(self.thefiles):
            tar = [0,0]
        else:
            tar = [0]
        rest = 0

        """The title of the whole figure gets the name of the Mesurement series"""
        if self.separate_borders_mode is True:
            fig.suptitle(self.all_vars[0]['Name']+  ' mode: '+mode+ 'separate_ntsr_int_border: %.3f' % self.ntsr_helper_border +comment, fontsize=9)
        else:
            fig.suptitle(self.all_vars[0]['Name']+' mode: '+mode+ 'ntsr_border: %.3f' % self.ntsr_border +comment, fontsize=9)
        for i in range(0, len(self.thefiles)):
            # print('tar ist '+ str(tar))
            x, y = np.meshgrid(np.arange(0, (self.all_vars[i]['Size']+1) * 0.39, 0.39),
                               np.arange(0, (self.all_vars[i]['Size']+1) * 0.39, 0.39))
            # z = np.reshape(np.array(all_integrals[i]), (all_vars[i]['Size'], all_vars[i]['Size']))
            if 'integrals' in mode:
                if 'nev_abs' in mode:
                    z = self.all_integrals_nev_abs[i]
                    bounds = [self.all_integrals_nev_abs.min(), self.all_integrals_nev_abs.max()]
                    axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm'+'\nnever absorbed power: ' + str(np.round(np.sum(z)*1E-6, 3)) + ' W', fontsize=8, style='italic')
                    # print('z pos: '+str(self.all_vars[i]['z_pos'])+'\nnever absorbed power: ' + str(np.round(np.sum(z)*1E-6, 3)) + ' W')
                elif 'remaining' in mode:
                    z = self.all_integrals_remaining[i]
                    bounds = [self.all_integrals_remaining.min(), self.all_integrals_remaining.max()]
                    axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm'+'\nremaining power: ' + str(np.round(np.sum(z)*1E-6, 3)) + ' W', fontsize=8, style='italic')
                    # print('z pos: '+str(self.all_vars[i]['z_pos'])+'\nremaining power: ' + str(np.round(np.sum(z)*1E-6, 3)) + ' W')
                else:
                    z = self.all_integrals[i]
                    bounds = [self.all_integrals.min(), self.all_integrals.max()]
                    axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm'+'\nreceived power: ' + str(np.round(np.sum(self.all_integrals[i])*1E-6, 3)) + ' W', fontsize=8, style='italic')
                    # print('z pos: '+str(self.all_vars[i]['z_pos'])+'\nreceived power: ' + str(np.round(np.sum(self.all_integrals[i])*1E-6, 3)) + ' W')
            elif 'ratios' in mode:
                z = self.all_ratios[i]
                bounds = [self.all_ratios.min(), self.all_ratios.max()]
                axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm', fontsize=8, style='italic')
            elif 'custom' in mode:
                z = data[i]
                bounds = [data.min(), data.max()]
                axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm', fontsize=8, style='italic')
            elif 'photons' in mode:
                if 'nev_abs' in mode:
                    z = self.all_photon_integrals_nev_abs[i]
                    bounds = [self.all_photon_integrals_nev_abs.min(), self.all_photon_integrals_nev_abs.max()]
                    axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm'+'\nnever absorbed photons: ' + str(np.round(np.sum(z)*1e-6, 3)) + ' \u03BCmol s$^{-1}$', fontsize=8, style='italic')
                elif 'remaining' in mode:
                    z = self.all_photon_integrals_remaining[i]
                    bounds = [self.all_photon_integrals_remaining.min(), self.all_photon_integrals_remaining.max()]
                    axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm'+'\nremaining photons: ' + str(np.round(np.sum(z)*1e-6, 3)) + ' \u03BCmol s$^{-1}$', fontsize=8, style='italic')
                else:
                    z = self.all_photon_integrals[i]
                    bounds = [self.all_photon_integrals.min(), self.all_photon_integrals.max()]
                    axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm'+'\nreceived photons: ' + str(np.round(np.sum(z)*1e-6, 3)) + ' \u03BCmol s$^{-1}$', fontsize=8, style='italic')

            if 'color' not in mode:
                if bar_to_max is False:
                    bounds = [z.min(), z.max()]
                    boundsstring = ''
                else:
                    boundsstring = '_boundstomax_'
                if 'contour_clabel' in mode:
                    if i == 0:
                        logging.info('contour_clabel mode detected')
                    # https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html
                    x_c, y_c = np.meshgrid(np.arange(0, (self.all_vars[i]['Size']) * 0.39, 0.39),
                               np.arange(0, (self.all_vars[i]['Size']) * 0.39, 0.39))
                    contour = axes[tuple(tar)].contour(x_c, y_c, z, 5, colors='black')
                    plt.clabel(contour, inline=True, fontsize=4, fmt='%i')
                    pcm = axes[tuple(tar)].pcolormesh(x, y, z, antialiased=False, shading='flat'
                                                      # , ec='face', alpha=0.5
                                                      )
                elif 'contourf' in mode:
                    if i == 0:
                        logging.info('contourf mode detected')
                    # https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html
                    x_c, y_c = np.meshgrid(np.arange(0, (self.all_vars[i]['Size']) * 0.39, 0.39),
                               np.arange(0, (self.all_vars[i]['Size']) * 0.39, 0.39))
                    pcm = axes[tuple(tar)].contourf(x_c, y_c, z, 20, cmap='viridis')
                elif 'contour' in mode:
                    if i == 0:
                        logging.info('contour mode detected')
                    # https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html
                    x_c, y_c = np.meshgrid(np.arange(0, (self.all_vars[i]['Size']) * 0.39, 0.39),
                               np.arange(0, (self.all_vars[i]['Size']) * 0.39, 0.39))
                    pcm = axes[tuple(tar)].contour(x_c, y_c, z, 20, cmap='viridis')
                else:
                    x_c, y_c = np.meshgrid(np.arange(0, (self.all_vars[i]['Size']) * 0.39 + 0.39, 0.39),
                                           np.arange(0, (self.all_vars[i]['Size']) * 0.39 + 0.39, 0.39))
                    pcm = axes[tuple(tar)].pcolormesh(x_c, y_c, z, antialiased=False, shading='flat',vmin=bounds[0], vmax=bounds[1])
            if 'color' in mode:
                axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm', fontsize=8, style='italic')
                z = self.all_colors[i]
                pcm = axes[tuple(tar)].imshow(z, origin='lower',extent=[x.min(), x.max(), y.min(), y.max()])
            yticks = np.linspace(0., 0.39 * (self.all_vars[i]['Size']), 7)
            yticks[1:-1] = np.round(np.linspace(0., 0.39 * (self.all_vars[i]['Size']), 7)[1:-1])
            xticks = np.linspace(0, 0.39 * (self.all_vars[i]['Size']), 7)
            xticks[1:-1] = np.round(np.linspace(0, 0.39 * (self.all_vars[i]['Size']), 7)[1:-1])
            axes[tuple(tar)].set_yticks(yticks)
            axes[tuple(tar)].set_yticks(xticks)
            divider = make_axes_locatable(axes[tuple(tar)])
            if 'color' not in mode:
                cax = divider.append_axes('right', size='5%', pad=0.05)
                cbar = fig.colorbar(pcm, cax=cax)
                cbar.ax.tick_params(labelsize=7, pad=0.1)
            if 'integrals' in mode:
                cbar.ax.set_ylabel("P / \u03BCW", fontsize=8, labelpad=0.2)
            elif 'ratios' in mode:
                cbar.ax.set_ylabel("peak ratio / 1", fontsize=8, labelpad=0.2)
            elif 'custom' in mode:
                cbar.ax.set_ylabel(cbar_title, fontsize=8, labelpad=0.2)
            elif 'photons' in mode:
                cbar.ax.set_ylabel('$q_\mathrm{n,p}$ / nmol s$^{-1}$', fontsize=8, labelpad=0.2)
            axes[tuple(tar)].tick_params(axis='x', labelsize=7, pad=2)
            axes[tuple(tar)].tick_params(axis='x', labelsize=7, pad=2)
            axes[tuple(tar)].tick_params(axis='y',labelsize=7, pad=0.1)
            axes[tuple(tar)].set_xlabel("x / cm", fontsize=8, labelpad=0.1)
            axes[tuple(tar)].set_ylabel("y / cm", fontsize=8, labelpad=0.1)
            axes[tuple(tar)].set_aspect('equal', adjustable='box')
            if self.separate_borders_mode is True:
                axes[tuple(tar)].text(0.1, y.max()+0.2, 'ntsr border: %.3f'%self.separate_ntsr[i, 0, 0], fontsize=4.5)
            if cols < len(self.thefiles):
                if tar[1] < cols-1:
                    tar[1] += 1
                elif tar[1] == cols-1:
                    tar[1] = 0
                    tar[0] +=1

                if i == len(self.thefiles)-1:
                    if len(self.thefiles) % cols == 0:
                        rest = 0
                    else:
                        rest = (len(self.thefiles)//cols+1)*cols-len(self.thefiles)
                    if rest != 0:
                        axes[tuple(tar)].remove()
                        rest -= 1
                while rest != 0:
                    if tar[1] < cols-1:
                        tar[1] += 1
                    elif tar[1] == cols-1:
                        tar[1] = 0
                        tar[0] += 1
                    axes[tuple(tar)].remove()
                    rest -= 1
            else:
                tar[0] += 1
                if len(self.thefiles) < cols and i == len(self.thefiles)-1:
                    axes[tuple(tar)].remove()
        fig.subplots_adjust(
                        left=0.125,  # the left side of the subplots of the figure
                        # right = 0.9,  # the right side of the subplots of the figure
                        #bottom = 0.1,  # the bottom of the subplots of the figure
                        top=0.95,     # the top of the subplots of the figure
                        wspace=0.4,  # the amount of width reserved for space between subplots,
                                      # expressed as a fraction of the average axis width
                        hspace=0.1)
        if self.thedir == '':
            self.thedir = Path(self.dir_path)
        fig.savefig(self.thedir/Path(self.all_vars[0]['Name'] +'Colormaps_%s'%mode+ comment+ boundsstring + str(len(self.thefiles)) + '.png')
                     , dpi=300
                     , bbox_inches='tight'
                    )
        plt.show()

    def plot_single(self, mode='integrals', data=False, info=False, comment = '', cbar_title='', bar_to_max=False):
        """
        Plots and saves one specific heatmap.
        -------------------
        Parameters:

            mode (string, default = 'integrals'):
                See docstring of method plot.

            data (2D numpy array):
                The array which should be plotted, should match paramaeter mode.

            comment (string, default = ''):
                Additional string to be displayed in the figure.

            info (dictionary):
                The dictionary containing the corresponding information to the given array. E.g. a subarray of all_vars.

            cbar_title (string, default = ''):
                Title of the colorbar for mode == 'custom'.

            bar_to_max (boolean, default is False):
                If True, all colorbars are scaled to the absolute maximum of the figure set.
        """
        fig, axes = plt.subplots(figsize=(3.4, 3))
        """The title of the whole figure becomes the name of the Mesurement series"""
        fig.suptitle(self.all_vars[0]['Name']+' mode: '+mode + 'ntsr_border: '+str(self.ntsr_border)+comment, fontsize=9)
        x, y = np.meshgrid(np.arange(0, (info['Size']+1) * 0.39, 0.39),
                           np.arange(0, (info['Size']+1) * 0.39, 0.39))
        if mode == 'integrals':
            z = data
            axes.set_title('z pos: '+str(info['z_pos'])+' cm'+'\nreceived power: ' + str(np.round(np.sum(z)*1E-6, 3)) + ' W', fontsize=8, style='italic')
        elif mode == 'ratios':
            z = data
            axes.set_title('z pos: '+str(info['z_pos'])+' cm', fontsize=8, style='italic')
        elif mode == 'custom':
            z = data
            axes.set_title('z pos: '+str(info['z_pos'])+' cm', fontsize=8, style='italic')
        if mode != 'color':
            pcm = axes.pcolormesh(x, y, z, antialiased=False, shading='flat')
        if mode == 'color':
            z = data
            pcm = axes.imshow(z, origin='lower',extent=[x.min(), x.max(), y.min(), y.max()])
        yticks = np.linspace(0., 0.39 * (info['Size']), 7)
        yticks[1:-1] = np.round(np.linspace(0., 0.39 * (info['Size']), 7)[1:-1])
        xticks = np.linspace(0, 0.39 * (info['Size']), 7)
        xticks[1:-1] = np.round(np.linspace(0, 0.39 * (info['Size']), 7)[1:-1])
        axes.set_yticks(yticks)
        axes.set_yticks(xticks)
        divider = make_axes_locatable(axes)
        if mode != 'color':
                bounds = [z.min(), z.max()]
                cax = divider.append_axes('right', size='5%', pad=0.05)
                if bar_to_max is True:
                    barticks = np.arange(bounds[0], bounds[1], round((bounds[1]-bounds[0])/6))
                    cbar = fig.colorbar(pcm, cax=cax, boundaries=bounds, ticks=barticks)
                else:
                    cbar = fig.colorbar(pcm, cax=cax)
                cbar.ax.tick_params(labelsize=7, pad=0.1)
        if mode == 'integrals':
            cbar.ax.set_ylabel("P / \u03BCW", fontsize=8, labelpad=0.2)
        elif mode == 'ratios':
            cbar.ax.set_ylabel("peak ratio / 1", fontsize=8, labelpad=0.2)
        elif mode == 'custom':
            cbar.ax.set_ylabel(cbar_title, fontsize=8, labelpad=0.2)
        axes.tick_params(axis='x', labelsize=7, pad=2)
        axes.tick_params(axis='y',labelsize=7, pad=0.1)
        axes.set_xlabel("x / cm", fontsize=8, labelpad=0.1)
        axes.set_ylabel("y / cm", fontsize=8, labelpad=0.1)
        axes.set_aspect('equal', adjustable='box')
        if self.thedir == '':
            self.thedir = Path(self.dir_path)
        fig.savefig(self.thedir/Path(info['Name'] +'Colormap_single_%sat%05.02f'%(mode,info['z_pos'])+ comment + str(len(self.thefiles)) + '.png')
                     , dpi=300
                     , bbox_inches='tight'
                    )
        plt.show()

    def to_minimize(self, value, percent, array):
        """The function to be minimized by the cutoff at percent function"""
        return abs(np.sum(np.where(array>=value, array, 0)) - np.sum(array)*percent/100)


    def cutoff_at_percent(self, percent, array):
        """
        Returns an array where the outer pixels are set to 0 to achieve a integral percentage value of each array
        measurement. Minimze results are stored in self.cutoff_results.
        -------------
        Parameters:

            percent (number):
                The percentage of the array to keep.

            array (array):
                The array to be treated.
        """
        self.cutoff_results = np.array([scipy.optimize.OptimizeResult]*len(array))
        cutoff_array = np.zeros(array.shape)
        for i in range(0, len(array)):
            self.cutoff_results[i] = minimize(self.to_minimize, np.max(array[i])*0.1, tol=1, args=(percent, array[i]), method= 'Powell')
            cutoff_array[i] = np.where(array[i] >= self.cutoff_results[i]['x'], array[i], 0)
        return cutoff_array

    def determine_angle(self, array, center = 'auto', comment=''):
        """
        Determines the emission angle depending on the input array. returns a plot and the calculated arrays. Also
        calculates the theretical orgin of the lightsource.
        ---------------
        Parameters:

            array (array):
                The array to be treated.

            center (tuple of two integers or 'auto', default is 'auto'):
                The center of the light source on the array. If center is ’auto’ , the position of the largest value in
                the array is used. Otherwise, pass a tuple of coordinates in pixels.

        ---------------
        Returns:

            angles (1D numpy array of 4 numbers):
                The opening angles in x and y direction.
        """
        if 'auto' in center:
            center = np.reshape(np.array(np.where(array[-1] == np.amax(array[-1]))), (2))
            logging.info('center found at ' + str(center)+ 'resp.'+ str(center*0.39))
        x1 = np.zeros(len(array))
        x2 = np.zeros(len(array))
        y1 = np.zeros(len(array))
        y2 = np.zeros(len(array))
        try:
            self.all_z[0] = self.all_z[0]
        except:
            self.all_z = np.zeros(len(self.all_vars))
            for num, i in np.ndenumerate(self.all_vars):
                self.all_z[num] =  self.all_vars[num]['z_pos']
        for i in range(0, len(array)):
            y1[i] = np.min(np.where(array[i, center[0], :] != 0))
            y2[i] = self.all_vars[-1]['Size']-np.min(np.where(np.flip(array[i, center[0], :]) != 0))
            x1[i] = np.min(np.where(array[i, :, center[1]] != 0))
            x2[i] = self.all_vars[-1]['Size']-np.min(np.where(np.flip(array[i, :, center[1]]) != 0))
        self.y1_model = np.polyfit(self.all_z, y1, 1)
        self.y2_model = np.polyfit(self.all_z, y2, 1)
        self.x1_model = np.polyfit(self.all_z, x1, 1)
        self.x2_model = np.polyfit(self.all_z, x2, 1)
        x, y, z = sp.symbols('x, y, z')
        eqy1 = sp.Eq(self.y1_model[0]*x + self.y1_model[1], y)
        eqy2 = sp.Eq(self.y2_model[0]*x + self.y2_model[1], y)
        self.ans_y = sp.solve((eqy1, eqy2), (x, y), set=True)[1].pop()
        logging.info('found origin (z) and center (x) from cut along x and y=%.2f: '%center[0]+ str(self.ans_y))
        eqx1 = sp.Eq(self.x1_model[0]*x + self.x1_model[1], y)
        eqx2 = sp.Eq(self.x2_model[0]*x + self.x2_model[1], y)
        self.ans_x = sp.solve((eqx1, eqx2), (x, y), set=True)[1].pop()
        logging.info('found origin (z) and center (y) from cut along y and y=%.2f: '%center[1]+ str(self.ans_x))
        angles = np.zeros(4)
        angles[0] = abs(math.degrees(math.sin(self.y1_model[0]*0.39)))
        angles[1] = abs(math.degrees(math.sin(self.y2_model[0]*0.39)))
        angles[2] = abs(math.degrees(math.sin(self.x1_model[0]*0.39)))
        angles[3] = abs(math.degrees(math.sin(self.x2_model[0]*0.39)))
        print(angles)
        print('mean: %f'%np.mean(angles))
        print('std: %f'%np.std(angles))
        fig, ax = plt.subplots(1, figsize=(3.3,2.8))
        ax.set_title(self.all_vars[0]['Name'], size=9)
        ax.set_xlabel('z / cm', fontsize=6, labelpad=0.1)
        ax.set_ylabel('x and y / cm', fontsize=6, labelpad=0.1)
        ax.tick_params(labelsize=6, pad=2)

        ax.plot(self.all_z, y1, 'ob')
        ax.plot(self.all_z, y2, 'ob')
        ax.plot(self.all_z, x1, 'or')
        ax.plot(self.all_z, x2, 'or')

        x_plot = np.append([self.ans_x[0]-2],self.all_z)
        y_plot = np.append([self.ans_y[0]-2],self.all_z)

        ax.plot(y_plot, self.y1_model[1]+self.y1_model[0]*y_plot, '-b')
        ax.plot(y_plot, self.y2_model[1]+self.y2_model[0]*y_plot, '-b')
        ax.plot(x_plot, self.x1_model[1]+self.x1_model[0]*x_plot, '-r')
        ax.plot(x_plot, self.x2_model[1]+self.x2_model[0]*x_plot, '-r')

        ax.plot(self.ans_y[0], self.ans_y[1],'*b')
        ax.plot(self.ans_x[0], self.ans_x[1],'*r')

        ax.text(0.99, 0.5, 'opening angle: %.2f deg pm %.2f'%(np.mean(angles), np.std(angles)), horizontalalignment='right',verticalalignment='center', transform=ax.transAxes, size=6)
        ax.text(0.02, 0.95, 'origin @ %.2f cm (z) \nand center @ %.2f (x) cut along x, y=%.2f\n'
                            'opening angle: %.2f deg pm %.2f'%(self.ans_y[0],self.ans_y[1],center[0],np.mean(angles[:2]), np.std(angles[:2])), horizontalalignment='left',verticalalignment='top', transform=ax.transAxes, wrap=True, size=5, c='b')
        ax.text(0.02, 0.05, 'origin @ %.2f cm (z) \nand center @ %.2f (y) cut along y, x=%.2f\n'
                           'opening angle: %.2f deg pm %.2f'%(self.ans_x[0],self.ans_x[1],center[1],np.mean(angles[2:]), np.std(angles[2:])), horizontalalignment='left',verticalalignment='bottom', transform=ax.transAxes, wrap=True, size=5, c='r')
        filename = self.thedir/Path(self.all_vars[0]['Name'] + '_find_angle_'+comment+'.png')
        plt.savefig(filename, dpi=600)
        plt.show()
        return angles

    def plot_slice_set(self, axis, borders=False, num_cuts=18, comment='', mode='integrals', info=False):
        """
        Plots a set of cuts through the chosen axis.
        ------------
        Parameters:

            axis (string):
               The axis on which the cut positions are set. Either 'x' or 'y'.

            borders (tuple of to integers, optional):
                The cuts are distributed according to a numpy linspace function. Here the respective start and stop can
                be set.

            num_cuts (integer, default = 18):
                The cuts are distributed according to a numpy linspace function. Here the respective number of cuts can
                be set.

            mode (string, default = 'integrals'):
                See docstring of method plot.

            comment (string, default = ''):
                Additional string to be displayed in the figure.

            info (dictionary):
                The dictionary containing the corresponding information to the given array. E.g. a subarray of all_vars.
        """
        if info is False:
            info = self.all_vars
        yticks = []
        for i in info:
            yticks.append(i['z_pos'])
        ygrids = np.zeros(len(yticks)+1)
        for i in range(0, len(yticks)+1):
            if i == 0:
                ygrids[i] = yticks[i]-(yticks[i+1]-yticks[i])/2
            elif i == len(yticks):
                ygrids[i] = yticks[i-1]+(yticks[i-1]-yticks[i-2])/2
            else:
                ygrids[i] = (yticks[i]+yticks[i-1])/2
        x, y = np.meshgrid(np.arange(0, (self.vars['Size']+1) * 0.39, 0.39), ygrids)
        """Here we set the size of the figure, when smaller than 6 we make two cols. when bigger, we make 3"""
        if num_cuts <= 6:
            cols = 2
        else:
            cols = 3
        fig, axes = plt.subplots(ncols=cols,
                                 nrows=math.ceil(num_cuts/cols),
                                 figsize=(cols*3.4, math.ceil(num_cuts/cols)*1),
                                 sharex=False,
                                 sharey=False)
        if cols < num_cuts:
            tar = [0, 0]
        else:
            tar = [0]
        rest = 0

        """The title of the whole figure gets the name of the Mesurement series"""
        fig.suptitle('$\mathit{sliceset~of:}$ \"'+self.all_vars[0]['Name']+'\" $\mathit{mode:}$ \"'+mode+ '\" $\mathit{ntsr\_border:}$ \"'+str(self.ntsr_border)+'\"'+comment, fontsize=9)
        if borders is False:
            cuts = np.linspace(0, info[0]['Size']*0.39, num_cuts+2)[1:-1]
        else:
            cuts = np.linspace(borders[0]*0.39, borders[1]*0.39, num_cuts)

        for cut in cuts:
            # print('cut: %.2f'%cut)
            # print('tar ist '+ str(tar))
            if 'integrals' in mode:
                if 'nev_abs' in mode:
                    if axis == 'y':
                        z = self.all_integrals_nev_abs[:, int(self.pos_to_num(cut)), :]
                    if axis == 'x':
                        z = self.all_integrals_nev_abs[:, :, int(self.pos_to_num(cut))]
                    axes[tuple(tar)].set_title('Slice of nev. abs. int. at %s = %.2f' % (axis, float(self.pos_to_num(cut)*0.39+0.39/2)), fontsize=8, style='italic')
                elif 'remaining' in mode:
                    if axis == 'y':
                        z = self.all_integrals_remaining[:, int(self.pos_to_num(cut)), :]
                    if axis == 'x':
                        z = self.all_integrals_remaining[:, :, int(self.pos_to_num(cut))]
                    axes[tuple(tar)].set_title('Slice of rem. int. at %s = %.2f' % (axis, float(self.pos_to_num(cut)*0.39+0.39/2)), fontsize=8, style='italic')
                else:
                    if axis == 'y':
                        z = self.all_integrals[:, int(self.pos_to_num(cut)), :]
                    if axis == 'x':
                        z = self.all_integrals[:, :, int(self.pos_to_num(cut))]
                    axes[tuple(tar)].set_title('Slice of integrals at %s = %.2f' % (axis, float(self.pos_to_num(cut)*0.39+0.39/2)), fontsize=8, style='italic')
            if 'photons' in mode:
                if 'nev_abs' in mode:
                    if axis == 'y':
                        z = self.all_photon_integrals_nev_abs[:, int(self.pos_to_num(cut)), :]
                    if axis == 'x':
                        z = self.all_photon_integrals_nev_abs[:, :, int(self.pos_to_num(cut))]
                    axes[tuple(tar)].set_title('Slice of nev. abs. phot. at %s = %.2f' % (axis, float(self.pos_to_num(cut)*0.39+0.39/2)), fontsize=8, style='italic')
                elif 'remaining' in mode:
                    if axis == 'y':
                        z = self.all_photon_integrals_remaining[:, int(self.pos_to_num(cut)), :]
                    if axis == 'x':
                        z = self.all_photon_integrals_remaining[:, :, int(self.pos_to_num(cut))]
                    axes[tuple(tar)].set_title('Slice of rem. phot. at %s = %.2f' % (axis, float(self.pos_to_num(cut)*0.39+0.39/2)), fontsize=8, style='italic')
                else:
                    if axis == 'y':
                        z = self.all_photon_integrals[:, int(self.pos_to_num(cut)), :]
                    if axis == 'x':
                        z = self.all_photon_integrals[:, :, int(self.pos_to_num(cut))]
                    axes[tuple(tar)].set_title('Slice of photons at %s = %.2f' % (axis, float(self.pos_to_num(cut)*0.39+0.39/2)), fontsize=8, style='italic')
            elif 'ratios' in mode:
                if axis == 'y':
                    z = self.all_ratios[:, int(self.pos_to_num(cut)), :]
                if axis == 'x':
                    z = self.all_ratios[:, :, int(self.pos_to_num(cut))]
                axes[tuple(tar)].set_title('Slice of ratios at %s = %.2f' % (axis, float(self.pos_to_num(cut)*0.39+0.39/2)), fontsize=8, style='italic')
            if 'color' not in mode:
                if 'contour_clabel' in mode:
                    if i == 0:
                        logging.info('contour_clabel mode detected')
                    # https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html
                    x_c, y_c = np.meshgrid(np.arange(0, (self.vars['Size']) * 0.39, 0.39), ygrids[:-1])
                    contour = axes[tuple(tar)].contour(x_c, y_c, z, 5, colors='black')
                    plt.clabel(contour, inline=True, fontsize=4, fmt='%i')
                    pcm = axes[tuple(tar)].pcolormesh(x, y, z, antialiased=False, shading='flat'
                                                      # , ec='face', alpha=0.5
                                                      )
                elif 'contourf' in mode:
                    if i == 0:
                        logging.info('contourf mode detected')
                    # https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html
                    x_c, y_c = np.meshgrid(np.arange(0, (self.vars['Size']) * 0.39, 0.39), ygrids[:-1])
                    pcm = axes[tuple(tar)].contourf(x_c, y_c, z, 20, cmap='viridis')
                elif 'contour' in mode:
                    if i == 0:
                        logging.info('contour mode detected')
                    # https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html
                    x_c, y_c = np.meshgrid(np.arange(0, (self.vars['Size']) * 0.39, 0.39), ygrids[:-1])
                    pcm = axes[tuple(tar)].contour(x_c, y_c, z, 20, cmap='viridis')
                else:
                    pcm = axes[tuple(tar)].pcolormesh(x, y, z, antialiased=False, shading='flat')
            if 'color' in mode:
                axes[tuple(tar)].set_title('Slice of colors at %s = %.2f' % (axis, float(self.pos_to_num(cut)*0.39+0.39/2)), fontsize=8, style='italic')
                if axis == 'y':
                    z = self.all_colors[:, self.pos_to_num(cut), :]
                if axis == 'x':
                    z = self.all_colors[:, :, self.pos_to_num(cut)]
                pcm = axes[tuple(tar)].imshow(z, origin='lower',extent=[x.min(), x.max(), y.min(), y.max()])
            xticks = np.linspace(0., 0.39 * (info[0]['Size']), 7)
            xticks[1:-1] = np.round(np.linspace(0., 0.39 * (info[0]['Size']), 7)[1:-1])
            # xticks[1:-1] = np.round(np.linspace(0, 0.39 * (info['Size']), 7)[1:-1])
            # axes[tuple(tar)].set_yticks(yticks)
            axes[tuple(tar)].set_xticks(xticks)
            divider = make_axes_locatable(axes[tuple(tar)])
            if 'color' not in mode:
                cax = divider.append_axes('right', size='5%', pad=0.05)
                cbar = fig.colorbar(pcm, cax=cax)
                cbar.ax.tick_params(labelsize=7, pad=0.1)
                """set the number of ticks"""
                cbar.ax.locator_params(nbins=3)
            if 'integrals' in mode:
                cbar.ax.set_ylabel("P / \u03BCW", fontsize=8, labelpad=0.2)
            elif 'ratios' in mode:
                cbar.ax.set_ylabel("peak ratio / 1", fontsize=8, labelpad=0.2)
            elif 'photons' in mode:
                cbar.ax.set_ylabel('$q_\mathrm{n,p}$ / nmol s$^{-1}$', fontsize=8, labelpad=0.2)
            # elif mode == 'custom':
            #     cbar.ax.set_ylabel(cbar_title, fontsize=8, labelpad=0.2)
            axes[tuple(tar)].tick_params(axis='x', labelsize=7, pad=2)
            axes[tuple(tar)].tick_params(axis='y',labelsize=7, pad=0.1)
            if axis == 'y':
                axes[tuple(tar)].set_xlabel("x / cm", fontsize=8, labelpad=0.1)
            if axis == 'x':
                axes[tuple(tar)].set_xlabel("y / cm", fontsize=8, labelpad=0.1)
            axes[tuple(tar)].set_ylabel("z / cm", fontsize=8, labelpad=0.1)
            axes[tuple(tar)].set_aspect('equal', adjustable='box')
            # axes[tuple(tar)].axis([x.min(), x.max(), y.min(), y.max()])
            # fig.colorbar(label="P / \u03BCW")
            # fig.tight_layout()


            if cols < num_cuts:
                if tar[1] < cols-1:
                    tar[1] += 1
                elif tar[1] == cols-1:
                    tar[1] = 0
                    tar[0] += 1
                if cut == cuts[-1]:
                    if num_cuts%cols == 0:
                        rest = 0
                    else:
                        rest = (num_cuts//cols+1)*cols-num_cuts
                    if rest != 0:
                        axes[tuple(tar)].remove()
                        rest -= 1
                while rest != 0:
                    if tar[1] < cols-1:
                        tar[1] += 1
                    elif tar[1] == cols-1:
                        tar[1] = 0
                        tar[0] +=1
                    axes[tuple(tar)].remove()
                    rest -= 1
            else:
                tar[0] += 1
                if num_cuts < cols and i == num_cuts-1:
                    axes[tuple(tar)].remove()
        fig.subplots_adjust(
                        left=0.125,  # the left side of the subplots of the figure
                        # right = 0.9,  # the right side of the subplots of the figure
                        #bottom = 0.1,  # the bottom of the subplots of the figure
                        top=0.95,     # the top of the subplots of the figure
                        wspace=0.4,  # the amount of width reserved for space between subplots,
                                      # expressed as a fraction of the average axis width
                        hspace=0.1)
        if self.thedir == '':
            self.thedir = Path(self.dir_path)
        fig.savefig(self.thedir/Path(info[0]['Name'] +'slice_set_%s_%s'%(mode, axis) + comment + str(len(self.thefiles)) + '.png')
                     , dpi=300
                     , bbox_inches='tight'
                    )
        plt.show()

    def plot_histograms(self, mode='ntsr', log=None, comment=''):
        """
        Plots and saves a set of histograms of the scans.
        -------------
        Parameters:

            mode (string, default is 'ntsr'):
                Available modes are:
                    'ntsr': Histograms of the NTSR values.
                    'peak_ratios': Histograms of the peak ratio values.
                    'integrals': Histograms of the received power per pixel.

            log (boolean, default is None):
                Defines if the y-axis of the histogram plot is in log scale. On default log is False for ’stds’ and True
                 for ’peak_ratios’ and ’integrals’.

            comment (string, optional):
                Will be displayed with the reactor name.
        """
        if len(self.thefiles) <= 6:
            cols = 2
        else:
            cols = 3
        fig, axes = plt.subplots(ncols=cols,
                                 nrows=math.ceil(len(self.thefiles)/cols),
                                 figsize=(cols*3.4, math.ceil(len(self.thefiles)/cols)*3),
                                 sharex=False,
                                 sharey=False)
        if cols < len(self.thefiles):
            tar = [0,0]
        else:
            tar = [0]
        rest = 0
        """The title of the whole figure gets the name of the Mesurement series"""
        if self.separate_borders_mode is True:
            fig.suptitle(self.all_vars[0]['Name'] + ' mode: ' + mode + ' separate_ntsr_int_border: %.3f' % self.ntsr_helper_border +comment, fontsize=9)
        else:
            fig.suptitle(self.all_vars[0]['Name'] + ' mode: '+mode + ' ntsr_border: %.3f' % self.ntsr_border +comment, fontsize=9)
        if mode == 'ntsr':
            data = self.all_ntsr
            if log is None:
                log = False
        elif mode == 'peak_ratios':
            data = self.all_ratios
            if log is None:
                log = True
        elif mode == 'integrals':
            data = self.all_integrals
            if log is None:
                log = True
        for i in range(0, len(self.thefiles)):
            axes[tuple(tar)].set_title('z pos: '+str(self.all_vars[i]['z_pos'])+' cm', fontsize=8, style='italic')
            axes[tuple(tar)].set_ylabel("frequency", fontsize=8, labelpad=0.1)
            if mode == 'ntsr':
                axes[tuple(tar)].set_ylabel("signal to noise / 1", fontsize=8, labelpad=0.1)
                n, bins, patches = axes[tuple(tar)].hist(self.all_ntsr[i].reshape(-1), 50, density=1, histtype='step', color=self.color_cycle_list[i], log=log)
                if self.separate_borders_mode is True:
                    axes[tuple(tar)].plot([self.separate_ntsr[i, 0, 0], self.separate_ntsr[i, 0, 0]], [0, n.max()/2], '-.'
                                    , label='%.3f'%self.separate_ntsr[i,0,0], alpha=0.5, color=self.color_cycle_list[i])
                    xmin, xmax = axes[tuple(tar)].get_xlim()
                    ymin, ymax = axes[tuple(tar)].get_ylim()
                    t_xpos = xmax+0.01*(xmax-xmin)
                    t_ypos = ymin+0.01*(ymax-ymin)
                    axes[tuple(tar)].text(t_xpos, t_ypos, 'ntsr_border = %.3f'%self.separate_ntsr[i,0,0], rotation='vertical', fontsize=6)
                if self.separate_borders_mode is False:
                    axes[tuple(tar)].plot([self.ntsr_border, self.ntsr_border], [0, n.max()], 'k-.'#, label='border'
                              )

            elif mode == 'peak_ratios':
                axes[tuple(tar)].set_ylabel("peak ratios / 1", fontsize=8, labelpad=0.1)
                n, bins, patches = axes[tuple(tar)].hist(self.all_ratios[i].reshape(-1), 50, density=1, histtype='step', color=self.color_cycle_list[i], log=log)

            elif mode == 'integrals':
                axes[tuple(tar)].set_ylabel("P / \u03BCW", fontsize=8, labelpad=0.1)
                n, bins, patches = axes[tuple(tar)].hist(self.all_integrals[i].reshape(-1), 50, density=1, histtype='step', color=self.color_cycle_list[i], log=log)

            if cols < len(self.thefiles):
                if tar[1] < cols-1:
                    tar[1] += 1
                elif tar[1] == cols-1:
                    tar[1] = 0
                    tar[0] +=1

                if i == len(self.thefiles)-1:
                    # rest = len(thefiles)%cols  # rest of empty axes, which are removed later
                    # rest = (tar[0]+1)**2-len(thefiles)
                    if len(self.thefiles)%cols == 0:
                        rest = 0
                    else:
                        rest = (len(self.thefiles)//cols+1)*cols-len(self.thefiles)
                        # print('rest ist '+str(rest))
                    if rest != 0:
                        axes[tuple(tar)].remove()
                        # print('removing '+str(tar))
                        rest -= 1
                while rest != 0:
                    if tar[1] < cols-1:
                        tar[1] += 1
                    elif tar[1] == cols-1:
                        tar[1] = 0
                        tar[0] +=1
                    # print('removing '+str(tar))
                    axes[tuple(tar)].remove()
                    rest -= 1
            else:
                tar[0] += 1
                # print('tar ist '+ str(tar))
                if len(self.thefiles) < cols and i == len(self.thefiles)-1:
                    axes[tuple(tar)].remove()
        fig.subplots_adjust(
                        left = 0.125,  # the left side of the subplots of the figure
                        # right = 0.9,  # the right side of the subplots of the figure
                        #bottom = 0.1,  # the bottom of the subplots of the figure
                        top = 0.95,     # the top of the subplots of the figure
                        wspace = 0.4,  # the amount of width reserved for space between subplots,
                                      # expressed as a fraction of the average axis width
                        hspace = 0.3)
        # fig.colorbar(label="P / \u03BCW")
        # fig.tight_layout()
        if self.thedir == '':
            self.thedir = Path(self.dir_path)
        fig.savefig(self.thedir/Path(self.all_vars[0]['Name'] +'_histograms_%s'%mode+ comment + str(len(self.thefiles)) + '.png')
                     , dpi=300
                     , bbox_inches='tight'
                    )
        plt.show()

    def plot_slice(self, axis, pos, comment='', mode='integrals', info=False, data=False, bar_to_max=False,
                   out=True, cbar_title='', for_str=False):
        """
        Plots and saves a cut through the chosen axis.
        --------------------
        Parameters:

            axis (string):
               The axis on which the cut positions are set. Either 'x' or 'y'.

            pos (number):
                Slice position in cm.

            mode (string, default = 'integrals'):
                See docstring of method plot.

            comment (string, default = ''):
                Additional string to be displayed in the figure.

            info (dictionary):
                The dictionary containing the corresponding information to the given array. E.g. a subarray of all_vars.

            data (array, default is False):
                Numpy array for mode == 'custom'.

            comment (string, default = ''):
                Additional string to be displayed in the figure.

            cbar_title (string, default = ''):
                Title of the colorbar for mode == 'custom'.

            bar_to_max (boolean, default is False):
                If True, all colorbars are scaled to the absolute maximum of the figure set.

            out (boolean, default is True):
                If True , the figure will be saved.

            for_str (string, optional):
                Defines the format of the cbar tick labels. E.g. ’:04.0f’. Handy for later animation to keep the size
                of the pictures similar.
        """
        if info is False:
            info = self.all_vars
        yticks = []
        for i in info:
            yticks.append(i['z_pos'])
        ygrids = np.zeros(len(yticks)+1)
        for i in range(0, len(yticks)+1):
            if i == 0:
                ygrids[i] = yticks[i]-(yticks[i+1]-yticks[i])/2
            elif i == len(yticks):
                ygrids[i] = yticks[i-1]+(yticks[i-1]-yticks[i-2])/2
            else:
                ygrids[i] = (yticks[i]+yticks[i-1])/2
        # print(yticks)
        # print(ygrids)
        x, y = np.meshgrid(np.arange(0, (self.vars['Size']+1) * 0.39, 0.39), ygrids)
        """Here we set the size of the figure, when smaller than 6 we make two cols. when bigger, we make 3"""
        fig, axes = plt.subplots(figsize=(3.4, 3))
        """The title of the whole figure gets the name of the Mesurement series"""
        fig.suptitle(self.all_vars[0]['Name']+' '+comment)
        # z = np.reshape(np.array(all_integrals[i]), (all_vars[i]['Size'], all_vars[i]['Size']))
        if mode == 'integrals':
            if axis == 'y':
                z = self.all_integrals[:, self.pos_to_num(pos), :]
            if axis == 'x':
                z = self.all_integrals[:, :, self.pos_to_num(pos)]
            bounds = [self.all_integrals.min(), self.all_integrals.max()]
            axes.set_title('Slice of integrals at %s = %.2f' % (axis, self.pos_to_num(pos)*0.39), fontsize=8, style='italic')
        elif mode == 'ratios':
            if axis == 'y':
                z = self.all_ratios[:, self.pos_to_num(pos), :]
            if axis == 'x':
                z = self.all_ratios[:, :, self.pos_to_num(pos)]
            bounds = [self.all_ratios.min(), self.all_ratios.max()]
            axes.set_title('Slice of ratios at %s = %.2f' % (axis, self.pos_to_num(pos)*0.39), fontsize=8, style='italic')
        elif mode == 'photons':
            if axis == 'y':
                z = self.all_photon_integrals[:, self.pos_to_num(pos), :]
            if axis == 'x':
                z = self.all_photon_integrals[:, :, self.pos_to_num(pos)]
            bounds = [self.all_photon_integrals.min(), self.all_photon_integrals.max()]
            axes.set_title('Slice of photons at %s = %.2f' % (axis, self.pos_to_num(pos)*0.39), fontsize=8, style='italic')
        elif mode == 'custom':
            if axis == 'y':
                z = data[:, self.pos_to_num(pos), :]
            if axis == 'x':
                z = data[ :, :, self.pos_to_num(pos)]
            bounds = [data.min(), data.max()]
            axes.set_title(' custom slice at %s = %.2f' % (axis, self.pos_to_num(pos)*0.39), fontsize=8, style='italic')
        if mode != 'color':
            if bar_to_max is False:
                bounds = [z.min(), z.max()]
            pcm = axes.pcolormesh(x, y, z, antialiased=False, shading='flat', vmin=bounds[0], vmax=bounds[1])
        if mode == 'color':
            if axis == 'y':
                z = self.all_colors[:, self.pos_to_num(pos), :]
            if axis == 'x':
                z = self.all_integrals[:, :, self.pos_to_num(pos)]
            pcm = axes.imshow(z, origin='lower',extent=[x.min(), x.max(), y.min(), y.max()])
        xticks = np.linspace(0., 0.39 * (info[0]['Size']), 7)
        xticks[1:-1] = np.round(np.linspace(0., 0.39 * (info[0]['Size']), 7)[1:-1])
        # xticks[1:-1] = np.round(np.linspace(0, 0.39 * (info['Size']), 7)[1:-1])
        # axes.set_yticks(yticks)
        axes.set_xticks(xticks)
        divider = make_axes_locatable(axes)
        if mode != 'color':
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cbar = fig.colorbar(pcm, cax=cax)
            cbar.ax.tick_params(labelsize=7, pad=0.1)
            if for_str is not False:
                cbar.ax.set_yticklabels([for_str.format(i) for i in cbar.get_ticks()]) # set ticks of your format
        if mode == 'integrals':
            cbar.ax.set_ylabel("P / \u03BCW", fontsize=8, labelpad=0.2)
        elif mode == 'ratios':
            cbar.ax.set_ylabel("peak ratio / 1", fontsize=8, labelpad=0.2)
        elif mode == 'photons':
            cbar.ax.set_ylabel("photons / nmol/s", fontsize=8, labelpad=0.2)
        elif mode == 'custom':
            cbar.ax.set_ylabel(cbar_title, fontsize=8, labelpad=0.2)
        axes.tick_params(axis='x', labelsize=7, pad=2)
        axes.tick_params(axis='y',labelsize=7, pad=0.1)
        if axis == 'y':
            axes.set_xlabel("x / cm", fontsize=8, labelpad=0.1)
        if axis == 'x':
            axes.set_xlabel("y / cm", fontsize=8, labelpad=0.1)
        axes.set_ylabel("z / cm", fontsize=8, labelpad=0.1)
        axes.set_aspect('equal', adjustable='box')
        if self.thedir == '':
            self.thedir = Path(self.dir_path)
        if out is True:
            fig.savefig(self.thedir/Path(info[0]['Name'] +'slice_single_%s_%s_at%04i'%(mode, axis, pos*100) + comment + str(len(self.thefiles)) + '.png')
                         , dpi=300
                         , bbox_inches='tight'
                        )
        plt.show()


    def residual(self, pars, x, LED, data):
        """Helper method to define the absorption of MeOr."""
        model = LED*(1-pars['stray_factor'])/(10**(skewed_gaussian(x, pars['amp'], pars['cen'], pars['wid'], pars['skew']))) + LED*pars['stray_factor']
        return model - data

    def out_all(self, slice=True, bar_to_max=False):
        """
        Outputs and saves all possible plot sets by the methods plot, plot_slice_set, SI_plot and plot_histograms.
        (62 files)
        --------------
        Parameters:

            slice (boolean, default is True):
                Defines if also all slicesets are plotted.

            bar_to_max (boolean, default is True):
                If True, all colorbars are scaled to the absolute maximum of the figure set.
        """
        self.SI_plot()
        plt.close()
        for i in ['ntsr', 'peak_ratios', 'integrals']:
            self.plot_histograms(mode=i)
            plt.close()
        for i in ['integrals_nev_abs', 'integrals_nev_abs_contour', 'integrals_nev_abs_contourf', 'integrals_nev_abs_contour_clabel',
            'integrals_remaining', 'integrals_remaining_contour', 'integrals_remaining_contourf', 'integrals_remaining_contour_clabel',
            'integrals', 'integrals_contour', 'integrals_contourf', 'integrals_contour_clabel',
            'photons_nev_abs', 'photons_nev_abs_contour', 'photons_nev_abs_contourf', 'photons_nev_abs_contour_clabel',
            'photons_remaining', 'photons_remaining_contour', 'photons_remaining_contourf', 'photons_remaining_contour_clabel',
            'photons', 'photons_contour', 'photons_contourf', 'photons_contour_clabel',
                  'ratios', 'ratios_contour', 'ratios_contourf', 'ratios_contour_clabel', 'color']:
            self.plot(mode=i, bar_to_max=bar_to_max)
            plt.close()
            if slice is True:
                self.plot_slice_set(axis='x', mode=i)
                plt.close()
                self.plot_slice_set(axis='y', mode=i)
                plt.close()


__author__ = "Maximilian Sender"
__copyright__ = 'Copyright 2021, Institute of Chemical Engineering, Ulm University'
__license__ = "GNU LGPL"
