"""
PyWKSpectra
===========

Wheeler-Kiladis Space-Time Spectral Analysis in Python.

Author: Alejandro Jaramillo-Moreno
Instituto de Ciencias de la Atmósfera y Cambio Climático
Universidad Nacional Autónoma de México (UNAM)
"""

from .io import load_wk_data
from .spectra import compute_wk_spectra
from .plotting import (
    plot_raw_spectra,
    plot_wk_spectra,
    plot_matsuno,
    matsuno_dispersion,
)

__version__ = "1.0"

__all__ = [
    "load_wk_data",
    "compute_wk_spectra",
    "plot_raw_spectra",
    "plot_wk_spectra",
    "plot_matsuno",
    "matsuno_dispersion",
]
