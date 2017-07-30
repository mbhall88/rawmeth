"""
This module is intended to handle plotting data for Samples and Fast5 files.
"""
from __future__ import print_function
import sys
import warnings
import itertools
from bokeh import palettes
from bokeh.io import output_file, save as bokeh_save
from bokeh.plotting import Figure
from bokeh.models.tickers import FixedTicker
from bokeh.models import FuncTickFormatter
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
        fig = Figure(
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
        col_map = palettes.all_palettes[col_key]
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
        maps = palettes.all_palettes.keys()
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

    @property
    def threshold(self):
        return self.arguments['threshold']

    @threshold.setter
    def threshold(self, threshold):
        self.arguments['threshold'] = threshold

    def create_plot(self, save_as=None):

        # create the figure object
        fig = self.figure()

        # if sample start loop through samples.
        for sample_idx, sample in enumerate(self.data):
            print('Plotting started for {}...      '.format(sample.name),
                  end='\r')
            sys.stdout.flush()
            x_data, y_data = self._get_plot_data(sample)

            fig.multi_line(x_data,
                           y_data,
                           line_width=self.line_width,
                           alpha=self.alpha,
                           color=self.colours[sample_idx],
                           legend=sample.name)
            print('Plotting finished for {}        '.format(sample.name))

        # else loop through fast5 files.

        fig = self._axis_formatting(fig)

    def _axis_formatting(self, fig):
        fig.xaxis.minor_tick_line_color = None
        fig.xaxis.axis_line_color = None
        fig.xaxis.bounds = (0.5, len(self.motif) - 0.5)
        fig.xaxis.major_tick_line_color = None

        ticks = [x + 0.5 for x in range(len(self.motif))]
        fig.xaxis.ticker = FixedTicker(ticks=ticks)

        # a dictionary to map the x position to it's motif nucleotide
        label_dict = {i + 0.5: label
                      for i, label in enumerate(self.motif)}

        # creates a function that will map the xaxis label to it's base.
        axis_formatter = FuncTickFormatter(
            code='var labels = {};'
                 'return labels[+tick];'.format(label_dict))
        fig.xaxis.formatter = axis_formatter

        return fig

    def _get_plot_data(self, reads):
        """Generates the x- and y-axis data for plotting.

        Args:
            reads (list[Fast5]): List of Fast5s to extract data from.

        Returns:
            (tuple[list[int|float], list[int|float]]): A tuple of the x and y
            data for the line plot.

        """
        x_data = []
        y_data = []

        for read_idx, read in enumerate(reads):
            motif_idxs = read.motif_indices(self.motif)

            for i in motif_idxs:
                signal_df = read.extract_motif_signal(i)

                if signal_df.empty:
                    continue

                if self.threshold:
                    signal_df = signal_df[signal_df[self.y_variable]
                        .between(self.threshold[0], self.threshold[1],
                                 inclusive=True)]

                x_data.append(_generate_line_plot_xs(signal_df['pos']))
                y_data.append(signal_df[self.y_variable])

        return x_data, y_data


def _generate_line_plot_xs(xs):
    """Generates x axis coordinates for producing a line plot of signals.

    As each event has a variable length, in order to plot a bunch of events
    together they need to be 'squashed' into the same x-axis domain.
    This is done by grouping all of the signals for a single event together,
    iterating through them and dividing their position in the event by the
    length of the event. This is effectively calculating the position in the
    event as a percentage of the way from the start to the end. NOTE: The
    position as a percentage is 0 indexed. i.e [0, 0] => [0.0, 0.5] not
    [0.5, 1.0]

    Examples:
        xs = [0, 0, 1, 1, 1, 1, 2, 2, 0, 0]
        _generate_line_plot_xs(xs)
        >>> [0.0, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 0.0, 0.5]

    Args:
        xs (list): A list of integers or floats.

    Returns:
        x_coords (list[float]): A list of the x-coordinates to use to
        'squash` the given events onto the same x-axis domain.

    """
    x_coords = []
    xs_generator = itertools.groupby(xs)
    for key, gen in xs_generator:
        event_len = sum(1 for _ in gen)
        pos_as_perc = [(float(pos) / event_len) + int(key)
                       for pos in range(event_len)]
        x_coords.extend(pos_as_perc)
    return x_coords
