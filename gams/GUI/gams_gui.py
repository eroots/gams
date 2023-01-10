from __future__ import annotations

import scipy.signal as signal
from skimage.filters import butterworth
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
from PyQt5.uic import loadUiType
from PyQt5 import QtWidgets, QtCore, QtGui
import numpy as np
from copy import deepcopy
from geoh5py.workspace import Workspace
from geoh5py.ui_json.input_file import InputFile
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import os
import sys
from gams.utils import padding
# from skimage.restoration import inpaint #/ Not enough use cases to justify a dependency on scikit-image


path = os.path.dirname(os.path.realpath(__file__))
Ui_MainWindow, QMainWindow = loadUiType(os.path.join(path, 'gams_gui.ui'))
window_title = 'Geoscience Analyst Magnetics Suite'


#####################################


class format_model_coords(object):
    # Custom class to return cursor coordinates and grid values
    def __init__(self, im, X, Y, x_label='y', y_label='y', data_units='nT', median=None, std=None):
        self.im = im
        self.x_label = x_label
        self.y_label = y_label
        self.X = X
        self.Y = Y
        if 'SI' in data_units:
            self.data_units = 'SI x 10^-3'
        else:
            self.data_units = data_units
        self.data_label = 'Value'
        self.median = median
        self.std = std

    def __call__(self, x, y):

        for ix, xx in enumerate(self.X):
            if xx > x:
                x_idx = ix - 1#min(ix, len(self.X) - 1) - 1
                break
            x_idx = len(self.X) - 2
        for iy, yy in enumerate(self.Y):
            if yy > y:
                y_idx = iy - 1#min(iy, len(self.Y) - 1) - 1
                break
            y_idx = len(self.Y) - 2
        vals = np.array(self.im.get_array())
        vals = np.reshape(vals, (len(self.X), len(self.Y)))[y_idx, x_idx]
        return '\t'.join(['{}: {:>4.4g} {}',
                          '{}: {:>4.4g} {}\n',
                          '{}: {:>4.4g} {}\n',
                          'Median: {:>4.4g}, STD: {:>4.4g}\n']).format(self.x_label, (y), '',
                                                                       self.y_label, (x), '',
                                                                       self.data_label, (vals), self.data_units,
                                                                       self.median, self.std)

class GamsViewer(QMainWindow, Ui_MainWindow):
    def __init__(self, input_file: InputFile = None):
        super(GamsViewer, self).__init__()
        # Dictionaries of the units / titles for plotted grids
        self.units = {'appsusc': r'SI $\mathregular{\times 10 ^{-3}}$',
                      'magtrans': 'nT',
                      'tiltAS0trans': 'radians',
                      'dxASA0': 'nT/m',
                      'dyASA0': 'nT/m',
                      'dhASA0': 'nT/m',
                      'dzASA0': 'nT/m',
                      'ASA0': 'nT',
                      'ASA1': 'nT/m',
                      'wavenumber': '1/km',
                      'original': 'nT',
                      'padded': 'nT',
                      'filtered': 'nT',
                      'tapered': 'nT',
                      'taper': '',
                      'fft': '',
                      'filtered_fft': '',
                      'extrapolated': 'nT'}
        self.titles = {'appsusc': 'Apparent Susceptiblity',
                       'magtrans': 'Magnetic Field at the Pole with Vertical Dip',
                       'tiltAS0trans': 'Tilt to Transform data to Pole and Vertical Dip',
                       'dxASA0': r'East Derivative of ASA${0}$',
                       'dyASA0': r'North Derivative of ASA${0}$',
                       'dhASA0': r'Horizontal Derivative of ASA${0}$',
                       'dzASA0': r'Vertical Derivative of ASA${0}$',
                       'ASA0': 'Zeroth-Order Analytic-Signal Amplitude',
                       'ASA1': 'First-Order Analytic-Signal Amplitude',
                       'wavenumber': 'Zeroth-Order Local Wavenumber',
                       'original': 'Magnetic Field',
                       'padded': 'Padded Mag. Field',
                       'filtered': 'Filtered Mag. Field',
                       'fft': 'Amplitude spectrum',
                       'filtered_fft': 'Filtered Amplitude Spectrum',
                       'tapered': 'Tapered Mag. Field',
                       'taper': 'Taper',
                       'extrapolated': 'Extrapolated Mag. Field'}
        self.filter_order = 3
        self.path = ''
        self.initialized = False
        self.setupUi(self)

        if input_file is None:
            self.load_workspace()
            self.load_grid()
        else:
            self.uijson_initialize(input_file)

        if self.initialized:
            self.setup_widgets()
            self.set_colourmap()
            self.calculate_grids(force=True)
            self.update_plots()

    def setup_widgets(self):
        # Initialize and connect widgets
        # Colour maps
        self.colour_group = QtWidgets.QActionGroup(self)
        self.colour_group.addAction(self.action_bwr)
        self.colour_group.addAction(self.action_jet)
        self.colour_group.addAction(self.action_hsv)
        self.colour_group.addAction(self.action_viridis)
        self.colour_group.addAction(self.action_grey)
        self.colour_group.setExclusive(True)
        self.colourmap = self.colour_group.checkedAction().text()
        self.colour_group.triggered.connect(self.set_colourmap)
        self.action_invert.triggered.connect(self.invert_cmap)
        self.action_colourbar.triggered.connect(self.update_plots)
        # Set up mpl window
        self.fig = Figure()
        self.subplots = [self.fig.add_subplot(111)]
        self.addmpl(self.fig)
        # Combo boxes
        self.combo_extrap.currentIndexChanged.connect(self.dummy_update_plots)
        self.combo_padding.currentIndexChanged.connect(self.dummy_update_plots)
        self.combo_taper.currentIndexChanged.connect(self.dummy_update_plots)
        self.combo_filter.currentIndexChanged.connect(self.dummy_update_plots)
        # Plot selections
        self.check_original.clicked.connect(self.update_plots)
        self.check_filtered.clicked.connect(self.update_plots)
        self.check_fft.clicked.connect(self.update_plots)
        self.check_filtered_fft.clicked.connect(self.update_plots)
        self.check_tapered.clicked.connect(self.update_plots)
        self.check_padded.clicked.connect(self.update_plots)
        self.check_extrapolated.clicked.connect(self.update_plots)
        self.check_taper.clicked.connect(self.update_plots)
        self.check_wavenumber.clicked.connect(self.update_plots)
        self.check_ASA0.clicked.connect(self.update_plots)
        self.check_ASA1.clicked.connect(self.update_plots)
        self.check_dzASA0.clicked.connect(self.update_plots)
        self.check_dxASA0.clicked.connect(self.update_plots)
        self.check_dyASA0.clicked.connect(self.update_plots)
        self.check_dhASA0.clicked.connect(self.update_plots)
        self.check_tiltAS0trans.clicked.connect(self.update_plots)
        self.check_magtrans.clicked.connect(self.update_plots)
        self.check_appsusc.clicked.connect(self.update_plots)
        # Calculation buttons
        self.button_recalculate.clicked.connect(self.recalculate)
        # Populate the default parameters
        self.extrapolation = self.combo_extrap.currentText().lower()
        self.padding = self.combo_padding.currentText().lower()
        self.taper = self.combo_taper.currentText().lower()
        # self.spin_extrap_param.setValue(int(self.nx / 2))
        # self.extrap_param = self.spin_extrap_param.value()
        self.taper_param = self.spin_taper_param.value()
        self.num_std = 2
        self.spin_low_cut.setMinimum(int(2*(self.dx + 1)))
        self.spin_low_cut.setMaximum(int(self.nx * self.dx))
        self.spin_low_cut.setValue(int(2*(self.dx + 1)))
        self.spin_high_cut.setMinimum(int(2*(self.dx + 2)))
        self.spin_high_cut.setMaximum(int(self.nx * self.dx))
        self.spin_high_cut.setValue(int(self.nx * self.dx))
        # self.spin_dc_buffer.setValue(0)
        # self.spin_dc_buffer.setMinimum(0)
        # self.spin_dc_buffer.setMaximum(int(self.nx * self.dx / 2))
        # self.spin_dc_buffer.setSingleStep(int(self.dx))
        self.low_cut = self.spin_low_cut.value()
        self.high_cut = self.spin_high_cut.value()
        # Spin box parameters
        self.spin_taper_param.valueChanged.connect(self.dummy_update_plots)
        self.spin_threshold.valueChanged.connect(self.dummy_update_plots)
        self.spin_upward_continuation.valueChanged.connect(self.dummy_update_plots)
        self.spin_field_intensity.valueChanged.connect(self.dummy_update_plots)
        self.spin_low_cut.valueChanged.connect(self.dummy_update_plots)
        self.spin_high_cut.valueChanged.connect(self.dummy_update_plots)
        # self.spin_dc_buffer.valueChanged.connect(self.dummy_update_plots)
        # self.spin_extrap_param.valueChanged.connect(self.dummy_update_plots)
        # Plot options check boxes
        self.action_link_axes.triggered.connect(self.update_plots)
        # self.check_tight_layout.clicked.connect(self.update_plots)
        # Open / Write menu actions
        self.action_workspace.triggered.connect(self.load_new_workspace)
        self.action_grid.triggered.connect(self.load_new_grid)
        self.action_write.triggered.connect(self.write_to_workspace)
        # Extra display options (menu)
        self.action_deviations.triggered.connect(self.set_num_deviation)
        self.action_deviation_cutoff.triggered.connect(self.update_plots)

    def load_workspace(self):
        # Opens a file dialog to load a new workspace (geoh5 file)
        while True:
            file, ret = QtWidgets.QFileDialog.getOpenFileName(self, 'Open Workspace', self.path,
                                                               'Geoscience Analyst (*.geoh5);; All Files (*)')
            if ret:
                self.workspace = Workspace(file)
                return True
            elif not (ret or self.initialized):
                ret = QtWidgets.QMessageBox.question(self, '', 'A geoh5 file is needed to run the GUI. Do you want to exit?',
                                                     QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                if ret == QtWidgets.QMessageBox.Yes:
                    sys.exit()
            else:
                return False

    def load_grid(self):
        # Opens a dialog to choose which grid to load from the current workspace.
        # Note that currently it loads all possible grids from the workspace, no guarantee they are magnetic grids
        objects = list(self.workspace.list_objects_name.values())
        file, ret = QtWidgets.QInputDialog.getItem(self, 'Select Grid',
                                                  'List of files', objects, 0, False)
        if ret:
            self.grid_name = file
            self.grid = self.workspace.get_entity(file)[0]
            self.grid_vals = {}
            try:
                vals = self.grid.get_data(file)[0].values
            except IndexError:
                fields = self.grid.get_data_list()
                print(fields)
                if 'Visual Parameters' in fields:
                    fields.remove('Visual Parameters')
                if len(fields) > 1:
                    field, ret = QtWidgets.QInputDialog.getItem(self, 'Select Field',
                                                      'List of fields', fields, 0, False)
                    vals = self.grid.get_data(field)[0].values
                else:
                    vals = self.grid.get_data(fields[0])[0].values
            vals[np.isnan(vals)] = 0
            self.centroids = self.grid.centroids
            self.dx = np.diff(self.centroids[:,0])[0]
            self.ny, self.nx = self.grid.shape
            self.grid_vals.update({'original': np.flipud(np.reshape(vals, [self.nx, self.ny]))})
            self.initialized = True

    def load_new_workspace(self):
        # Load new workspace and grid, and recalculate + plot
        ret = self.load_workspace()
        if ret:
            self.load_grid()
            self.calculate_grids(force=True)
            self.update_plots()

    def load_new_grid(self):
        # Load a new grid from current workspace, recalculate + plot
        self.load_grid()
        self.calculate_grids(force=True)
        self.update_plots()

    def apply_extrapolation(self, force=False):
        # Extrapolate grid to fill holes and square off edges.

        # Don't recalculate if nothing has changed, unless force is True
        if (self.extrapolation != self.combo_extrap.currentText().lower()) or force:
        # or (self.extrap_param != self.spin_extrap_param.value()) or force:
            self.extrapolation = self.combo_extrap.currentText().lower()
            # self.extrap_param = self.spin_extrap_param.value()
        else:
            return False
        if self.extrapolation == 'median fill':
            grid_vals = deepcopy(self.grid_vals['original'])
            grid_vals[self.zeros_idx] = np.median(self.grid_vals['original'])
            self.grid_vals['extrapolated'] = grid_vals
        elif self.extrapolation == 'zero fill':
            grid_vals = deepcopy(self.grid_vals['original'])
            grid_vals[self.zeros_idx] = 0
            self.grid_vals['extrapolated'] = grid_vals
        elif self.extrapolation == 'infill':
            # Infill can be very time / resource intensive - keep the zeros to a minimum
            # I've actually disabled the infill option. All the relevant code can be kept
            # but its not really worth the required import for the small number of use cases
            # Ideally the grids being worked on should be in good shape beforehand, which could include inpainting.
            if sum(self.zeros_idx) > 100: # Parameter could be changed
                QtWidgets.QMessageBox.question(self, '', 'Number of zeros is >100; Inpainting could take a while. Proceed?',
                                                     QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                if ret == QtWidgets.QMessageBox.No:
                    self.extrapolation = 'linear'
                    self.combo_extra.setIndex(0)
            self.grid_vals['extrapolated'] = inpaint.inpaint_biharmonic(grid_vals, self.zeros_idx)
        elif self.extrapolation in ('linear', 'nearest', 'cubic'):
            if np.any(self.grid_vals['original'] == 0):
                extrap_multiple = 2
                pad_x = int(self.nx * extrap_multiple / 3)
                pad_y = int(self.ny * extrap_multiple / 3)
                extrap_vals = np.pad(self.grid_vals['original'], 
                                 [[pad_x, pad_x], [pad_y, pad_y]], mode='reflect')
                X, Y = np.meshgrid(np.linspace(0, self.ny+2*pad_y, self.ny+2*pad_y),
                                   np.linspace(0, self.nx+2*pad_x, self.nx+2*pad_x))
                idx = extrap_vals != 0
                Xp, Yp = X[idx].flatten(), Y[idx].flatten()
                vals = extrap_vals[idx].flatten()
                new_vals = griddata((Xp, Yp), vals, (X.flatten(), Y.flatten()), method=self.extrapolation)
                new_vals = np.reshape(new_vals, [self.nx+2*pad_x, self.ny+2*pad_y])
                self.grid_vals['extrapolated'] = new_vals[pad_x:pad_x+self.nx, pad_y:pad_y+self.ny]
            else:
                self.grid_vals['extrapolated'] = deepcopy(self.grid_vals['original'])
        return True

    def apply_padding(self, force=False):
        # Applies the selected padding method. Largely a wrapper for numpy padding methods
        # Don't recalculate if nothing has changed, unless force is True (or if going to or from PSD tapering)
        if (self.padding != self.combo_padding.currentText().lower()) or \
           (self.combo_taper.currentText().lower() == 'psd') or (self.taper == 'psd') \
            or force:
            self.padding = self.combo_padding.currentText().lower()
        else:
            return False
        if self.padding == 'zeros':
            self.grid_vals.update({'padded': padding.zeros(self.grid_vals['extrapolated'],
                                        [[self.pad_left, self.pad_right], [self.pad_bot, self.pad_top]])})
            # QtWidgets.QMessageBox.question(self, '', '{}'.format(self.grid_vals['padded'].shape))
        elif self.padding == 'wrap':
            self.grid_vals.update({'padded': padding.wrap(self.grid_vals['extrapolated'],
                                        [[self.pad_left, self.pad_right], [self.pad_bot, self.pad_top]])})
        elif self.padding == 'reflect':
            self.grid_vals.update({'padded': padding.reflect(self.grid_vals['extrapolated'],
                                        [[self.pad_left, self.pad_right], [self.pad_bot, self.pad_top]])})
        elif self.padding == 'reflect (inverse)':
            self.grid_vals.update({'padded': padding.reflect(self.grid_vals['extrapolated'],
                                        [[self.pad_left, self.pad_right], [self.pad_bot, self.pad_top]], inverse=True)})
        elif self.padding == 'derivative mirror':
            pad_vals = padding.derivative_mirror(self.grid_vals['extrapolated'],
                                                 'lr', [self.pad_left, self.pad_right, 0, 0])
            pad_vals = padding.derivative_mirror(pad_vals,
                                                 'ud', [0, 0, self.pad_bot, self.pad_top])
            self.grid_vals.update({'padded': pad_vals})
        return True

    def calculate_taper(self, force=False):
        # Calculate the taper, store it to apply later
        # Don't recalculate anything unless force is True
        if (self.padding != self.combo_taper.currentText().lower()) \
         or (self.taper_param != self.spin_taper_param.value()) or force:
            self.taper = self.combo_taper.currentText().lower()
            self.taper_param = self.spin_taper_param.value()
        else:
            return False
        if self.taper == 'psd':
            taper_vals = np.ones((self.padded_size, self.padded_size))
            # Periodic-Smooth Decomposition
            # Based on code from https://github.com/jacobkimmel/ps_decomp/blob/master/psd.py
            # self.grid_vals.update({'tapered': padding.psd(self.grid_vals['filtered'], padding=None)})
            self.grid_vals.update({'tapered': padding.psd(self.grid_vals['padded'],
                                                         padding=None)})
        elif self.taper == 'blur':
            # Blur doesn't actually use a taper - set it to ones
            taper_vals = np.ones((self.padded_size, self.padded_size))
        else:
            if self.taper.lower() == 'tukey':
                # Calculate the area to keep flat - should be a bounding rectangle around the orignal grid
                beta_y = 1 - (self.nx / (self.nx + self.pad_left + self.pad_right))
                beta_x = 1 - (self.ny / (self.ny + self.pad_bot + self.pad_top))
                t1 = np.vstack([signal.windows.get_window((self.taper, beta_x), self.padded_size)])
                t2 = np.vstack([signal.windows.get_window((self.taper, beta_y), self.padded_size)]).T
            elif self.taper.lower() == 'kaiser':
                beta_x, beta_y = self.taper_param, self.taper_param
                t1 = np.vstack([signal.windows.get_window((self.taper, beta_x), self.padded_size)])
                t2 = np.vstack([signal.windows.get_window((self.taper, beta_y), self.padded_size)]).T
            else:
                t1 = np.vstack([signal.windows.get_window(self.taper, self.padded_size)])
                t2 = np.vstack([signal.windows.get_window(self.taper, self.padded_size)]).T
                                
            t1 = np.tile(t1, [self.padded_size, 1])
            t2 = np.tile(t2, [1, self.padded_size])
            taper_vals = (t1 * t2)
        self.grid_vals.update({'taper': taper_vals})
        return True

    def apply_taper(self):
        # Apply the calculated taper. 
        # Unless its set to blur - then calculate and apply the blur
        if self.taper == 'blur':
            binary = deepcopy(self.grid_vals['padded'])
            binary[abs(binary)>0] = 1
            # binary[self.zeros_idx] = 0
            med = deepcopy(self.grid_vals['padded'])
            med[self.grid_vals['padded']==0] = np.median(self.grid_vals['padded'][self.grid_vals['padded']!=0])
            blur = gaussian_filter(med, sigma=self.spin_taper_param.value())
            mask = gaussian_filter(binary, sigma=self.spin_taper_param.value())
            blur[mask>=1] = self.grid_vals['padded'][mask>=1]
            blur[blur == 0] = np.median(self.grid_vals['original'][np.where(self.grid_vals['original'] != 0)])
            self.grid_vals['taper'] = mask
            self.grid_vals['tapered'] = deepcopy(blur)
        elif self.taper == 'psd':
            return
            # self.grid_vals.update({'tapered': self.grid_vals['padded'] * self.grid_vals['taper']})
        else:
            self.grid_vals.update({'tapered': self.grid_vals['padded'] * self.grid_vals['taper']})

    def apply_filter(self, force=False):
        # Apply the chosen wavelength filter (if any)
        # dc_buffer = self.dx * self.nx - self.spin_dc_buffer.value()
        if self.spin_low_cut.value() >= self.spin_high_cut.value():
            self.spin_low_cut.setValue(self.low_cut)
            self.spin_high_cut.setValue(self.high_cut)
        else:
            self.low_cut = self.spin_low_cut.value()
            self.high_cut = self.spin_high_cut.value()
        self.grid_vals['filtered'] = deepcopy(self.grid_vals['tapered'])
        self.grid_vals['fft'] = abs(np.fft.fftshift(np.fft.fft2(self.grid_vals['filtered'])))
        self.grid_vals['filtered_fft'] = deepcopy(self.grid_vals['fft'])
        if self.combo_filter.currentText().lower() == 'none':
            if force:
                return True
            else:
                return False
        #######################
        # Naive FFT-based filtering approach

        # if self.combo_filter.currentText().lower() == 'lowpass':
        #     b, a = signal.butter(self.filter_order, [1/self.low_cut], fs=1/self.dx, btype='low')
        # elif self.combo_filter.currentText().lower() == 'highpass':
        #     b, a = signal.butter(self.filter_order, [1/self.high_cut], fs=1/self.dx, btype='high')
        #     b2, a2 = signal.butter(self.filter_order, [1/(dc_buffer)], fs=1/self.dx, btype='low')
        # elif self.combo_filter.currentText().lower() == 'bandpass':
        #     b, a = signal.butter(self.filter_order, [1/self.high_cut, 1/self.low_cut], fs=1/self.dx, btype='band')
        #     b2, a2 = signal.butter(self.filter_order, [1/(dc_buffer)], fs=1/self.dx, btype='low')
        # w, h = signal.freqz(b, a, fs=1/self.dx, worN=self.padded_size, whole=True)
        # w1 = np.tile(np.abs(h), [self.padded_size, 1])
        # w2 = np.tile(np.abs(h), [self.padded_size, 1]).T
        # filt = np.fft.fftshift(w1 * w2)
        # self.grid_vals['fft'] = np.fft.fftshift(np.fft.fft2(self.grid_vals['tapered']))
        # self.grid_vals['filtered_fft'] = self.grid_vals['fft'] * filt
        # # Add back in the DC component
        # if self.combo_filter.currentText().lower() in ['highpass', 'bandpass'] and self.spin_dc_buffer.value() != 0:
        #     w, h = signal.freqz(b2, a2, fs=1/self.dx, worN=self.padded_size, whole=True)
        #     # h = np.flip(h)
        #     h = 1 - h
        #     w1 = np.tile(np.abs(h), [self.padded_size, 1])
        #     w2 = np.tile(np.abs(h), [self.padded_size, 1]).T
        #     filt = np.fft.fftshift(w1 * w2)
        #     dc_component = self.grid_vals['fft'] * filt
        #     self.grid_vals['filtered_fft'] += dc_component
        # self.grid_vals['filtered'] = np.real(np.fft.ifft2(np.fft.ifftshift(self.grid_vals['filtered_fft'])))
        # self.grid_vals['filtered_fft'] = np.abs(self.grid_vals['filtered_fft'])
        # self.grid_vals['fft'] = np.abs(self.grid_vals['fft'])
        ##########################
        # Filtering using skimage.filters.butterworth
        # Also flip high/low to make things clearer for the user
        # 'High'
        low_cut = self.low_cut / self.dx
        high_cut = self.high_cut / self.dx
        if self.combo_filter.currentText().lower() == 'lowpass':
            self.grid_vals['filtered'] = butterworth(self.grid_vals['tapered'],
                                                     cutoff_frequency_ratio=1/low_cut,
                                                     high_pass=False,
                                                     order=self.filter_order)
        elif self.combo_filter.currentText().lower() == 'highpass':
            self.grid_vals['filtered'] = butterworth(self.grid_vals['tapered'],
                                                     cutoff_frequency_ratio=1/high_cut,
                                                     high_pass=True,
                                                     order=self.filter_order)
        elif self.combo_filter.currentText().lower() == 'bandpass':
            low_pass = butterworth(self.grid_vals['tapered'],
                                   cutoff_frequency_ratio=1/low_cut,
                                   high_pass=False,
                                   order=self.filter_order)
            high_pass = butterworth(self.grid_vals['tapered'],
                                    cutoff_frequency_ratio=1/high_cut,
                                    high_pass=True,
                                    order=self.filter_order)
            self.grid_vals['filtered'] = self.grid_vals['tapered'] - low_pass - high_pass
        self.grid_vals['fft'] = np.abs(np.fft.fftshift(np.fft.fft2(self.grid_vals['tapered'])))
        self.grid_vals['filtered_fft'] = np.abs(np.fft.fftshift(np.fft.fft2(self.grid_vals['filtered'])))
        
        # self.grid_vals['fft'] = np.abs(self.grid_vals['fft'])

    def recalculate(self):
        # Recalculate grids and plot them
        self.calculate_grids()
        self.update_plots()

    def calculate_grids(self, force=False):
        # Pre-process data and calculate the transformed data from it
        self.setWindowTitle(window_title + ' - Calculating extrapolation...')
        force = self.apply_extrapolation(force)
        self.setWindowTitle(window_title + ' - Calculating padding...')
        force = self.apply_padding(force)
        # force = self.apply_filter(force)
        self.setWindowTitle(window_title + ' - Calculating taper...')
        self.calculate_taper(force)
        self.apply_taper()
        force = self.apply_filter(force)
        self.setWindowTitle(window_title + ' - Calculating fields...')
        # Calculate all those fields
        next_pow2 = self.padded_size
        tapered_vals = self.grid_vals['filtered']
        opVDFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        opHXFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        opHYFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        opHilXFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        opHilYFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        magFD = np.fft.fft2(self.grid_vals['filtered'])
        freq = 2.*np.pi*np.fft.fftfreq(next_pow2, d=self.dx)
        # block_freq is a tiled frequency matrix so we can avoid loops
        block_freq = np.tile(freq, [next_pow2, 1])
        opVDFD = np.sqrt(block_freq*block_freq + block_freq.T*block_freq.T)
        ######################
        opHXFD = np.tile(-1j*freq, [next_pow2, 1])
        opHYFD = np.tile(-1j*freq, [next_pow2, 1]).T
        # opHYFD = np.tile(-1j*freq, [next_pow2, 1])
        # opHXFD = np.tile(-1j*freq, [next_pow2, 1]).T
        ######################
        # Preemptively zeroing the terms causes problems with where the values end up.
        idx = opVDFD == 0
        opVDFD[idx] = 1
        ######################
        # opHilYFD[idx] = -1j*block_freq[idx] / opVDFD[idx]
        # opHilXFD[idx] = -1j*block_freq[idx].T / opVDFD[idx]
        # opHilXFD[idx] = -1j*block_freq[idx] / opVDFD[idx]
        # opHilYFD[idx] = -1j*block_freq[idx].T / opVDFD[idx]
        opHilXFD = -1j*block_freq / opVDFD
        opHilYFD = -1j*block_freq.T / opVDFD
        ######################
        opHilYFD[opVDFD == 0.0] = 0.0
        opHilXFD[opVDFD == 0.0] = 0.0
        opVDFD[idx] = 0
        # print(magFD.shape)
        # print(opVDFD.shape)
        VDFD=magFD*opVDFD
        HXFD=magFD*opHXFD
        HYFD=magFD*opHYFD
        HilXFD=magFD*opHilXFD
        HilYFD=magFD*opHilYFD
        VD=np.real(np.fft.ifft2(VDFD))
        HX=np.real(np.fft.ifft2(HXFD))
        HY=np.real(np.fft.ifft2(HYFD))
        HilX=np.real(np.fft.ifft2(HilXFD))
        HilY=np.real(np.fft.ifft2(HilYFD))
        ASA1=np.sqrt(VD*VD+HX*HX+HY*HY)
        # Upward continue the field
        magUPFD = magFD * np.exp(-opVDFD * self.spin_upward_continuation.value())
        # Calculate hilbert transforms of this field
        HilXFDUP = magUPFD * opHilXFD
        HilYFDUP = magUPFD * opHilYFD
        # Transform back to space domain
        magUP = np.real(np.fft.ifft2(magUPFD))
        HilXUP = np.real(np.fft.ifft2(HilXFDUP))
        HilYUP = np.real(np.fft.ifft2(HilYFDUP))
        ASA0 = np.sqrt(tapered_vals*tapered_vals + HilX*HilX + HilY*HilY)
        ASA0UP = np.sqrt(magUP*magUP + HilXUP*HilXUP + HilYUP*HilYUP)
        dzASA0 = (ASA0UP - ASA0) / self.spin_upward_continuation.value()
        # Calculate the horizontal derivative of the ASA0 -- assume that it is zero at the outer edges
        dxASA0=np.zeros((next_pow2,next_pow2))
        dyASA0=np.zeros((next_pow2,next_pow2))
        ######################
        # dxASA0[1:-1, :] = (ASA0[2:, :] - ASA0[:-2, :]) / (self.dx*2)
        # dyASA0[:, 1:-1] = (ASA0[:, 2:] - ASA0[:, :-2]) / (self.dx*2)
        dyASA0[1:-1, :] = (ASA0[2:, :] - ASA0[:-2, :]) / (self.dx*2)
        dxASA0[:, 1:-1] = (ASA0[:, 2:] - ASA0[:, :-2]) / (self.dx*2)
        ######################
        # Calculate total horizontal derivative
        dhASA0 = np.sqrt(dxASA0*dxASA0 + dyASA0*dyASA0)
        tiltAS0trans = np.zeros((next_pow2,next_pow2))
        softnum = np.zeros((next_pow2,next_pow2))
        softdenom = np.zeros((next_pow2,next_pow2))
        thr = self.spin_threshold.value()
        coef = deepcopy(dzASA0)
        idx = abs(coef) > thr
        softnum[idx] = np.sign(coef[idx]) * (abs(coef[idx]) - thr)
        coef = abs(deepcopy(dhASA0))
        idx = abs(coef) > thr
        softdenom[idx] = coef[idx] - thr
        tiltAS0trans = np.arctan2(softnum, softdenom)
        magtrans = ASA0 * (np.sin(-tiltAS0trans))
        appsusc = 1000 * magtrans / (self.spin_field_intensity.value()*2.*np.pi)

        # Properly store some of these. Others could be stored, but these are the ones that are plottable.
        # Any new plotables would have to be added here.
        self.grid_vals.update({'appsusc': appsusc})
        self.grid_vals.update({'magtrans': magtrans})
        self.grid_vals.update({'tiltAS0trans': tiltAS0trans})
        self.grid_vals.update({'dhASA0': dhASA0})
        self.grid_vals.update({'dzASA0': dzASA0})
        self.grid_vals.update({'dxASA0': dxASA0})
        self.grid_vals.update({'dyASA0': dyASA0})
        self.grid_vals.update({'ASA0': ASA0})
        self.grid_vals.update({'ASA1': ASA1})
        self.grid_vals.update({'wavenumber': -(dzASA0 / ASA0) * 1000})
        # self.grid_vals.update({'VD': VD})
        # Trim the calculated grids back down to the original size
        # For debugging it could be interesting to view the full transformed fields, but not needed for general use.
        for key, val in self.grid_vals.items():
            if key not in ('original', 'tapered', 'extrapolated', 'padded', 'taper', 'fft', 'filtered_fft'):
                trimmed_vals = val[self.pad_left:self.pad_left + self.nx, self.pad_bot:self.pad_bot+self.ny]
                trimmed_vals[self.zeros_idx] = np.nan
                self.grid_vals.update({key: trimmed_vals})
        self.setWindowTitle(window_title)

    def calculate_grids_loops(self):
        # This is the original(ish) version of the code. Loops are very slow.
        self.apply_extrapolation()
        self.apply_padding()
        self.calculate_taper()
        self.apply_taper()
        # Calculate all those fields
        next_pow2 = self.padded_size
        tapered_vals = self.grid_vals['tapered']
        opVDFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        opHXFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        opHYFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        opHilXFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        opHilYFD=np.zeros((next_pow2,next_pow2),dtype=complex)
        magFD = np.fft.fft2(self.grid_vals['tapered'])
        freq = 2.*np.pi*np.fft.fftfreq(next_pow2, d=self.dx)
        for ii in range(0,next_pow2):
            for jj in range(0,next_pow2):
                opVDFD[ii,jj] = np.sqrt(freq[ii]*freq[ii]+freq[jj]*freq[jj])
                opHXFD[ii,jj] = -1j*freq[jj]
                opHYFD[ii,jj] = -1j*freq[ii]
                # opHXFD[ii,jj]=np.complex128(0.0-1j*freq[jj])
                # opHYFD[ii,jj]=np.complex128(0.0-1j*freq[ii])
                if(opVDFD[ii,jj]==0.0):
                    opHilXFD[ii,jj]=0.0
                    opHilYFD[ii,jj]=0.0
                else:
                    opHilYFD[ii,jj] = -1j*freq[ii]/opVDFD[ii,jj]
                    opHilXFD[ii,jj] = -1j*freq[jj]/opVDFD[ii,jj]
                    opHilYFD[ii,jj]=np.complex128(0.0-1j*freq[ii])/opVDFD[ii,jj]
                    opHilXFD[ii,jj]=np.complex128(0.0-1j*freq[jj])/opVDFD[ii,jj]

        VDFD=magFD*opVDFD
        HXFD=magFD*opHXFD
        HYFD=magFD*opHYFD
        HilXFD=magFD*opHilXFD
        HilYFD=magFD*opHilYFD
        VD=np.real(np.fft.ifft2(VDFD))
        HX=np.real(np.fft.ifft2(HXFD))
        HY=np.real(np.fft.ifft2(HYFD))
        HilX=np.real(np.fft.ifft2(HilXFD))
        HilY=np.real(np.fft.ifft2(HilYFD))
        ASA1=np.sqrt(VD*VD+HX*HX+HY*HY)
        
        #upward continue the field 10 metres
        # Is the 10 m general, or should be changable by user?
        # Make it changeable
        magUPFD=magFD*np.exp(-opVDFD*10.)
        #calculate hilbert transforms of this field
        HilXFDUP=magUPFD*opHilXFD
        HilYFDUP=magUPFD*opHilYFD
        #transform back to space domain
        magUP=np.real(np.fft.ifft2(magUPFD))
        HilXUP=np.real(np.fft.ifft2(HilXFDUP))
        HilYUP=np.real(np.fft.ifft2(HilYFDUP))
        ASA0=np.sqrt(tapered_vals*tapered_vals + HilX*HilX + HilY*HilY)
        ASA0UP=np.sqrt(magUP*magUP + HilXUP*HilXUP + HilYUP*HilYUP)
        dzASA0=(ASA0UP - ASA0)/10.
        #calculate the horizontal derivative of the ASA0 -- assume that it is zero at the outer edges
        dxASA0=np.zeros((next_pow2,next_pow2))
        dyASA0=np.zeros((next_pow2,next_pow2))
           
        for ii in range(1,next_pow2 - 1):
            for jj in range(1,next_pow2 - 1):
                dyASA0[jj,ii]=(ASA0[jj+1,ii]-ASA0[jj-1,ii])/(self.dx*2)
                dxASA0[ii,jj]=(ASA0[ii,jj+1]-ASA0[ii,jj-1])/(self.dx*2)
        # check the ASA0 and derivatives are calculated corectly
        #calculate total horizontal derivative
        dhASA0=np.sqrt(dxASA0*dxASA0 + dyASA0*dyASA0)
        tiltAS0trans=np.zeros((next_pow2,next_pow2))
        thr = self.spin_threshold.value()  #threshold value not required for this dataset.
        for ii in range (0, next_pow2, 1):
            for jj in range (0, next_pow2, 1):
            # threshold small values of the gradient
                coef=dzASA0[ii,jj]
                if abs(coef)>thr:
                    softnum=np.sign(coef)*(abs(coef)-thr)
                else:
                    softnum=0.0
                coef=abs(dhASA0[ii,jj])
                if abs(coef)>thr:
                    softdenom=(coef)-thr
                else:
                    softdenom=0.0
                tiltAS0trans[ii,jj]=np.arctan2(softnum,softdenom)
        magtrans= ASA0*(np.sin(-tiltAS0trans))
        appsusc= magtrans/(50000*2.*np.pi)

        # Properly store some of these. Others could be stored, but these are the ones that are plottable.
        self.grid_vals.update({'appsusc': appsusc})
        self.grid_vals.update({'magtrans': magtrans})
        self.grid_vals.update({'tiltAS0trans': tiltAS0trans})
        self.grid_vals.update({'dxASA0': dxASA0})
        self.grid_vals.update({'dyASA0': dyASA0})
        self.grid_vals.update({'dhASA0': dhASA0})
        self.grid_vals.update({'dzASA0': dzASA0})
        self.grid_vals.update({'ASA0': ASA0})
        self.grid_vals.update({'ASA1': ASA1})
        self.grid_vals.update({'wavenumber': -(dzASA0 / ASA0) * 1000})

    def invert_cmap(self):
        # Flip the selected colour map.
        if self.action_invert.isChecked():
            self.colourmap += '_r'
        else:
            self.colourmap = self.colourmap[:-2]
        self.update_plots()

    def set_colourmap(self):
        # Set the colour map based on checked menu option
        new_cmap = self.colour_group.checkedAction().text()
        if new_cmap != self.colourmap:
            self.colourmap = new_cmap
            if self.action_invert.isChecked():
                self.colourmap += '_r'
            self.update_plots()

    def set_num_deviation(self):
        val, ok_pressed = QtWidgets.QInputDialog.getDouble(self,
                                                           '# Standard Deviations',
                                                           'Value: ',
                                                           self.num_std, 0.1, 100, 2)
        if ok_pressed:
            self.num_std = val
            self.update_plots()

    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)  # Make a canvas
        self.mplvl.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(canvas=self.canvas,
                                         parent=self.mplwindow, coordinates=True)
        self.toolbar.setFixedHeight(48)

        self.canvas.draw()
        self.mplvl.addWidget(self.toolbar)

    def gen_tiling(self):
        # Generate an optimal plot grid for the number of active plots
        num_plots = len(self.active_plots)
        if num_plots < 3:
            tiling = [num_plots, 1]
        else:
            s1 = np.floor(np.sqrt(num_plots))
            s2 = np.ceil(num_plots / s1)
            tiling = [int(s1), int(s2)]
        return tiling

    def init_axes(self, tiling):
        # Initialize the axes based on a generating tiling
        self.axes = []
        for ii in range(tiling[0] * tiling[1]):
            if self.action_link_axes.isChecked() and ii > 0:
                self.axes.append(self.fig.add_subplot(tiling[0], tiling[1], ii + 1,
                                 sharex=self.axes[0], sharey=self.axes[0]))
            else:
                self.axes.append(self.fig.add_subplot(tiling[0], tiling[1], ii + 1))

    def dummy_update_plots(self, skip=False):
        # Dummy method to handle some of the widget signals - will also update plots if auto-update is on
        if self.combo_taper.currentText().lower() == 'blur':
            # self.padding = 'zeros'
            self.combo_padding.setCurrentIndex(0)
            self.combo_extrap.setCurrentIndex(1)
        if self.combo_taper.currentText().lower() in ('kaiser', 'blur'):
            self.spin_taper_param.setEnabled(True)
        else:
            self.spin_taper_param.setEnabled(False)
        if self.auto_update:
            self.calculate_grids()
            self.update_plots()

    def update_plots(self):
        # Clear and redraw the plots
        # Could be faster if only changed plots were re-drawn, but not a big problem right now.
        tiling = self.gen_tiling()
        self.fig.clear()
        self.init_axes(tiling)
        to_plot = self.active_plots
        self.plots = []
        for ii, grid in enumerate(to_plot):
            grid_vals = self.grid_vals[grid].flatten()
            median = np.nanmedian(grid_vals)
            std = np.nanstd(grid_vals)
            if self.action_deviation_cutoff.isChecked():
                vmin, vmax = (median - std*self.num_std), (median + std*self.num_std)
            else:
                vmin, vmax = (np.min(grid_vals), np.max(grid_vals))
            self.plots.append(self.axes[ii].imshow(self.grid_vals[grid], cmap=self.cmap,
                                                   vmin=vmin, vmax=vmax,
                                                   interpolation='none'))
            self.axes[ii].set_title(self.titles[grid])
            self.axes[ii].format_coord = format_model_coords(self.plots[ii],
                                                             X=list(range(self.grid_vals[grid].shape[0])),
                                                             Y=list(range(self.grid_vals[grid].shape[1])),
                                                             x_label='Row', y_label='Col.',
                                                             data_units=self.units[grid],
                                                             median=median, std=std)
            # Only give labels to the left- and bottom-most axes
            if ii % tiling[1] == 0:
                self.axes[ii].set_ylabel('Grid Row')
            if ii + 1 > tiling[1] * (tiling[0] - 1):
                self.axes[ii].set_xlabel('Grid Column')
            if self.action_colourbar.isChecked():
                cb = self.fig.colorbar(ax=self.axes[ii], mappable=self.plots[ii])
                cb.set_label(self.units[grid], rotation=270,
                                               labelpad=10,
                                               fontsize=12)
        # Remove any dangling axes
        while ii < tiling[0] * tiling[1] - 1:
            self.axes[ii+1].remove()
            ii += 1
        # if self.check_tight_layout.checkState():
        # self.fig.set_tight_layout(self.check_tight_layout.checkState() > 0)
        # print(self.check_tight_layout.checkState() > 0)
        self.canvas.draw()

    def uijson_initialize(self, input_file: InputFile):
        """
        Initialize the UI from ui.json file

        Parameters
        ----------
        input_file: InputFile object with parameters

        Returns
        -------

        """
        self.workspace = input_file.workspace
        self.grid_name = input_file.data["grid"].name
        self.grid = input_file.data["grid"]
        self.grid_vals = {}
        vals = input_file.data["data"].values
        vals[np.isnan(vals)] = 0
        self.centroids = self.grid.centroids
        self.dx = np.diff(self.centroids[:, 0])[0]
        self.ny, self.nx = self.grid.shape
        self.grid_vals.update({'original': np.flipud(np.reshape(vals, [self.nx, self.ny]))})
        self.initialized = True

    def write_to_workspace(self):
        # Write active plot grids to the current workspace.
        # If those grids already exist, by default they will be written with sequential numbering.
        # Note that unless you have GA Pro, you will have to reload the workspace in GA to access the new grids
        ret = QtWidgets.QMessageBox.question(self, 'Write to Workspace', 'OK to write all currently plotted fields to workspace?',
                                             QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if ret == QtWidgets.QMessageBox.Yes:
            names = self.active_plots
            fields = [self.grid_vals[n] for n in names]
            names = [self.titles[n] for n in names if n not in ('tapered', 'padded')]
            to_add = {n: {'values': np.flipud(f)} for n, f in zip(names, fields)}
            self.grid.add_data(to_add)
    
    @property
    def active_plots(self):
        # Get the currently checked plots
        plots = []
        if self.check_original.checkState():
            plots.append('original')
        if self.check_extrapolated.checkState():
            plots.append('extrapolated')
        if self.check_padded.checkState():
            plots.append('padded')
        if self.check_tapered.checkState():
            plots.append('tapered')
        if self.check_taper.checkState():
            plots.append('taper')
        if self.check_filtered.checkState():
            plots.append('filtered')
        if self.check_fft.checkState():
            plots.append('fft')
        if self.check_filtered_fft.checkState():
            plots.append('filtered_fft')
        if self.check_ASA0.checkState():
            plots.append('ASA0')
        if self.check_ASA1.checkState():
            plots.append('ASA1')
        if self.check_dzASA0.checkState():
            plots.append('dzASA0')
        if self.check_dxASA0.checkState():
            plots.append('dxASA0')
        if self.check_dyASA0.checkState():
            plots.append('dyASA0')
        if self.check_dhASA0.checkState():
            plots.append('dhASA0')
        if self.check_wavenumber.checkState():
            plots.append('wavenumber')
        if self.check_tiltAS0trans.checkState():
            plots.append('tiltAS0trans')
        if self.check_magtrans.checkState():
            plots.append('magtrans')
        if self.check_appsusc.checkState():
            plots.append('appsusc')
        return plots

    @property
    def auto_update(self):
        return self.check_auto_update.checkState()

    @property
    def cmap(self):
        return cm.get_cmap(self.colourmap)

    @property
    def padded_size(self):
        # The size of the padded grid. Should be a square unless Periodic-Smooth Decomp. is being used.
        if self.combo_taper.currentText().lower() == 'psd':
            return max(self.nx, self.ny)
        else:
            return int(2*(max(self.nx, self.ny)))

    @property
    def pad_left(self):
        return int(np.floor((self.padded_size - self.nx) / 2))

    @property
    def pad_right(self):
        return int(np.ceil((self.padded_size - self.nx) / 2))

    @property
    def pad_bot(self):
        return int(np.floor((self.padded_size - self.ny) / 2))

    @property
    def pad_top(self):
        return int(np.ceil((self.padded_size - self.ny) / 2))

    @property
    def ylim(self):
        return [self.pad_left + self.nx, self.pad_left]

    @property
    def xlim(self):
        return [self.pad_bot, self.pad_bot + self.ny]

    @property
    def zeros_idx(self):
        return self.grid_vals['original'] == 0

def main():
    app = QtWidgets.QApplication(sys.argv)

    ifile = None
    if len(sys.argv) > 1:
        ifile = InputFile.read_ui_json(sys.argv[1])

    viewer = GamsViewer(input_file=ifile)
    viewer.setWindowTitle(window_title)
    viewer.show()
    ret = app.exec_()
    sys.exit(ret)
    # viewer.disconnect_mpl_events()


if __name__ == '__main__':
    main()