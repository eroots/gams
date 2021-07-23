import numpy as np
from scipy.interpolate import griddata
# from skimage.restoration import inpaint


def linear(grid_vals):
	return grid_vals


def cubic(grid_vals):
	return grid_vals


def infill(grid_vals):
	return grid_vals
	# zeros_idx = grid_vals == 0
	# inpaint.inpaint_biharmonic(grid_vals, zeros_idx)
