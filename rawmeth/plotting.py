"""
This module is intended to handle plotting data for Samples and Fast5 files.
"""
from __future__ import print_function
import warnings
import bokeh
# from bokeh import palettes,
# from bokeh.io import output_file, save as bokeh_save
# from bokeh.plotting import Figure
# from bokeh.models.tickers import FixedTicker
# from bokeh.models import FuncTickFormatter
from rawmeth import fast5


DEFAULT_ARGS = {
    'plot_height': 960,
    'plot_width': 500,
    'title': None,
    'y_axis_label': None,
    'x_axis_label': None,
    'tools': 'pan, box_zoom, reset, save, wheel_zoom',
    'legend': True,
    'alpha': None,
    'linewidth': None,
    'colour_map': 'Set1',
    'threshold': None
}


class Plot(object):

    """This is a base class that is not intended for direct use.

    This class handles parsing of the arguments for plotting and setting up the
    figure environment.

    """

    def __init__(self, **kwargs):
        """Initialises a base class to handle argument parsing and setting up
        Figure dimesions etc. for plotting.

        Args:
            plot against each other. Can also just pass a single instance of
            either of these two classes.
            **kwargs (dict): Any arguments the user wishes to specify that
            differ from the defaults.
            TODO: add docs for defaults

        """
        self.valid_args = DEFAULT_ARGS.keys()[:]
        self.arguments = self._parse_kwargs(**kwargs)

    def __repr__(self):
        return self.figure()

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
    def height(self):
        """Returns the height used for the plot in pixels."""
        return self.arguments['plot_height']

    @height.setter
    def height(self, height):
        """Sets the height for the plot in pixels.

        Args:
            height (int): Height in pixels for the plot.

        """
        self.arguments['plot_height'] = height

    @property
    def width(self):
        """Returns the width used for the plot in pixels."""
        return self.arguments['plot_width']

    @width.setter
    def width(self, width):
        """Sets the width for the plot in pixels.

        Args:
            width (int): Width in pixels for the plot.

        """
        self.arguments['plot_width'] = width

    @property
    def title(self):
        """Returns the title for the plot."""
        return self.arguments['title']

    @title.setter
    def title(self, title):
        """Sets the title for the plot.

        Args:
            title (str): The desired title.
        """
        self.arguments['title'] = title

    @property
    def ylabel(self):
        """Returns the label that will be used for the y-axis."""
        return self.arguments['y_axis_label']

    @ylabel.setter
    def ylabel(self, label):
        """Sets the label to be used on the y-axis.

        Args:
            label (str): The desired y-axis label.
        """
        self.arguments['y_axis_label'] = label

    @property
    def xlabel(self):
        """Returns the label that will be used for the x-axis."""
        return self.arguments['x_axis_label']

    @xlabel.setter
    def xlabel(self, label):
        """Sets the label to be used on the y-axis.

        Args:
            label (str): The desired x-axis label.
        """
        self.arguments['x_axis_label'] = label

    @staticmethod
    def default_arguments():
        """Returns the dictionary of default arguments for plotting.

        Returns:
            DEFAULT_ARGS (dict): The default arguments dictionary.
        """
        return DEFAULT_ARGS

    def figure(self):
        """Creates and returns the base Bokeh Figure object that will be used
        by all the other plots in this module.

        Returns:
            (bokeh.plotting.Figure): The figure object that plots will be built
             off.

        """
        fig = bokeh.plotting.Figure(
            plot_height=self.height,
            plot_width=self.width,
            title=self.title,
            x_axis_label=self.xlabel,
            y_axis_label=self.ylabel,
            tools=self.arguments['tools']
        )

        if self.arguments['legend']:
            fig.legend.location = 'top_right'
            fig.legend.click_policy = 'hide'  # clicking group will hide it
            fig.legend.label_text_font = 'roboto'
            fig.legend.background_fill_alpha = 0

        return fig


class LinePlot(Plot):

    def __init__(self, data, motif, y_variable, **kwargs):
        self.motif = fast5.Motif(motif)
        self.data = data
        self.y_variable = y_variable
        super(LinePlot, self).__init__(**kwargs)

    @property
    def colours(self):
        """The colour map being used for plotting as a list of string
        representations of the hexidecimal RGB.

        Returns:
            (list[str]): a list of string representations of the
            hexidecimal RGB.

        """
        col_key = self.arguments['colour_map']
        col_map = bokeh.palettes.all_palettes[col_key]
        colours = col_map[col_map.keys()[-1]]
        return colours

    @colours.setter
    def colours(self, col_map):
        """Sets the colour map for the plot.

        Args:
            col_map (str): The name of the colour map scheme to use. Check
            http://bokeh.pydata.org/en/latest/docs/reference/palettes.html for
            the possible colour maps.

        """
        maps = bokeh.palettes.all_palettes.keys()
        if col_map not in maps:
            raise ValueError('Colour map must be one of: {}'.format(maps))
        self.arguments['colour_map'] = col_map

    @property
    def line_width(self):
        return self.arguments['linewidth']

    @line_width.setter
    def line_width(self, linewidth):
        self.arguments['linewidth'] = linewidth

    @property
    def alpha(self):
        return self.arguments['alpha']

    @alpha.setter
    def alpha(self, level):
        if not 0 <= level <= 1.:
            raise ValueError('Alpha must be between 0 and 1.')
        self.arguments['alpha'] = level
