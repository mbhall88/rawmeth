"""
This module is intended to handle plotting data for Samples and Fast5 files.
"""
from __future__ import print_function
import warnings
from bokeh import palettes
from bokeh.io import output_file, save as bokeh_save
from bokeh.plotting import Figure
from bokeh.models.tickers import FixedTicker
from bokeh.models import FuncTickFormatter
from rawmeth.fast5 import *


DEFAULT_ARGS = {
    'yaxis': 'signal',
    'alpha': 0.05,
    'linewidth': 1.,
    'colour_map': 'Set1',
    'save_as': None,
    'figsize': (960, 500),
    'theshold': (0, 2000),
    'legend': True
}


class Plot(object):

    """This is a base class that is not intended for direct use.

    This class handles parsing of the arguments for plotting and setting up the
    figure environment.

    """

    def __init__(self, motif, data, **kwargs):
        """Initialises a base class to handle argument parsing and setting up
        Figure dimesions etc. for plotting.

        Args:
            motif (Motif | str): Motif to plot.
            data (list[Sample | Fast5]): A list of samples or fast5 files to
            plot against each other. Can also just pass a single instance of
            either of these two classes.
            **kwargs (dict): Any arguments the user wishes to specify that
            differ (or not) from the defaults.

        """
        self.valid_args = DEFAULT_ARGS.keys()[:]
        self.arguments = self._parse_kwargs(**kwargs)
        self.motif = Motif(motif)

        # handles the case where user just passes a single object
        if isinstance(data, list):
            self.data = data
        else:
            self.data = [data]

    def _parse_kwargs(self, **kwargs):
        """Returns a dictionary to be used for plotting based on the
        arguments the user desires.

        Args:
            **kwargs (dict): A dictionary of parameters for plotting that the
            user wishes to use.

        Returns:
            (dict): A dictionary of the parameters to use to generate the plot.

        """
        args = DEFAULT_ARGS[:]

        for key, value in kwargs.items():
            if key not in self.valid_args:
                warnings.warn('{} is not a valid plotting argument. This '
                              'argument will be ignored.\nValid '
                              'arguments: {}'.format(key, self.valid_args))
                continue
            args[key] = value

        return args

    @property
    def colours(self):
        """The colour map being used for plotting as a list of string
        representations of the hexidecimal RGB.

        Returns:
            (list[str]): a list of string representations of the
            hexidecimal RGB.

        """
        col_key = self.arguments['colour_map']
        col_map = palettes.all_palettes[col_key]
        self.colours = col_map[col_map.keys()[-1]]
        return self.colours

    @colours.setter
    def colours(self, col_map):
        """Sets the colour map for the plot.

        Args:
            col_map (list[str]): a list of string representations of the
            hexidecimal RGB.

        """
        self.colours = col_map

    @property
    def height(self):
        """Returns the height used for the plot in pixels."""
        return self.arguments['figsize'][1]

    @height.setter
    def height(self, height):
        """Sets the height for the plot in pixels.

        Args:
            height (int): Height in pixels for the plot.

        """
        self.arguments['figsize'] = (self.width, height)

    @property
    def width(self):
        """Returns the width used for the plot in pixels."""
        return self.arguments['figsize'][0]

    @width.setter
    def width(self, width):
        """Sets the width for the plot in pixels.

        Args:
            width (int): Width in pixels for the plot.

        """
        self.arguments['figsize'] = (width, self.height)

    @property
    def title(self):
        """Returns the title for the plot."""
        self.title = 'Nanopore signal across {} motif'.format(self.motif)
        return self.title

    @title.setter
    def title(self, title):
        """Sets the title for the plot.

        Args:
            title (str): The desired title.
        """
        self.title = title

    @property
    def ylabel(self):
        """Returns the label that will be used for the y-axis."""
        yaxis = self.arguments['yaxis']
        if yaxis == 'signal':
            self.ylabel = 'Raw Signal (pA)'
        else:
            self.ylabel = 'Normalised Signal'
        return self.ylabel

    @ylabel.setter
    def ylabel(self, label):
        """Sets the label to be used on the y-axis.

        Args:
            label (str): The desired y-axis label.
        """
        self.ylabel = label

    @staticmethod
    def default_arguments():
        """Returns the dictionary of default arguments for plotting.

        Returns:
            DEFAULT_ARGS (dict): The default arguments dictionary.
        """
        return DEFAULT_ARGS
