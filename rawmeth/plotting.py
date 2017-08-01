"""
This module is intended to handle plotting data for Samples and Fast5 files.
"""
from __future__ import print_function
import sys
import warnings
import itertools
from copy import deepcopy
from bokeh import palettes
from bokeh.io import output_file, save as bokeh_save
from bokeh.plotting import Figure
from bokeh.models.tickers import FixedTicker
from bokeh.models import FuncTickFormatter
from rawmeth import fast5


DEFAULT_ARGS = {
    'plot_height': 600,
    'plot_width': 960,
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
        self.valid_args = deepcopy(DEFAULT_ARGS)
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
        args = deepcopy(DEFAULT_ARGS)

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

        return fig


class LinePlot(Plot):
    """Produces a line plot of the raw signal events related to the given
     motif.



    """
    def __init__(self, data, motif, y_variable, **kwargs):
        """Initialise a line plot for plotting the signal of fast5 files within
        a sample over a given motif.

        Args:
            data (list[Sample]): A list of Samples to plot.
            motif (Motif): The motif to plot signal for.
            y_variable (str): Variable to plot on the y-axis.
            **kwargs (dict): Arguments to format the plot. These include:
                plot_height (int): Height of plot in pixels.
                plot_width (int): Width of plot in pixels.
                title (str): Titlt for the plot.
                y_axis_label (str): Label for the y-axis
                x_axis_label (str): Label for the x-axis
                tools (str): Tools to add to the plot HTML. See
                http://bokeh.pydata.org/en/latest/docs/user_guide/tools.html
                for options.
                legend (bool): Whether or not to add a legend to the plot.
                alpha (float|int): The alpha (opacity) of the line. Must be
                between 0 and 1.
                linewidth (float|int): The width of the lines.
                colour_map (str): The colour scheme to use. See
                http://bokeh.pydata.org/en/latest/docs/reference/palettes.html
                for a list of available sets.
                threshold (tuple(int, int)): Limit the y-axis to the given
                threshold. The first element is the lower bound, second is the
                upper.

        """
        self.motif = fast5.Motif(motif)

        if not isinstance(data, list):
            self.data = [data]
        else:
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
        """Returns the current line width for the plot."""
        return self.arguments['linewidth']

    @line_width.setter
    def line_width(self, linewidth):
        """Sets the line width for the plot.

        Args:
            linewidth (int|float): Line width value for the lines.

        """
        self.arguments['linewidth'] = linewidth

    @property
    def alpha(self):
        """Returns the current alpha level for the plot."""
        return self.arguments['alpha']

    @alpha.setter
    def alpha(self, level):
        """Sets the alpha (opacity) level for the plot.

        Args:
            level (int|float): Alpha level between 0 and 1.

        """
        if not 0 <= level <= 1.:
            raise ValueError('Alpha must be between 0 and 1.')
        self.arguments['alpha'] = level

    @property
    def threshold(self):
        """Returns the current threshold for the plot."""
        return self.arguments['threshold']

    @threshold.setter
    def threshold(self, threshold):
        """Sets the current threshold level for the plot. This will filter the
        signals those whose values fall within this threshold.

        Args:
            threshold (tuple[int, int]): The lower (first element) and upper
            (second element) bounds for the threshold.

        """
        self.arguments['threshold'] = threshold

    def create_plot(self, save_as=None):
        """Creates the figure object that can be used to view the plot in a
        notebook environment, or to save to a file. If save_as is provided, the
        plot will be saved to a HTML file. Within this file you can interact
        with the plot (zoom, pan etc.) and save to PNG from the HTML view.

        Args:
            save_as (str): The file path to save the plot as. No file extension
            is required as .html will be added to the end of whatever argument
            is given.

        Returns:
            fig (bokeh.plotting.Figure): A figure object that can be viewed in
            an interactive environment such as a notebook. Or additional
            settings on the plot can be manipulated (using the Bokeh
            documentation).

        """
        # create the figure object
        fig = self.figure()

        # loop through the sample and generate the data to create the lines for
        # each fast5 file in the sample.
        for sample_idx, sample in enumerate(self.data):
            print('Plotting started for {}...      '.format(sample.name),
                  end='\r')
            sys.stdout.flush()

            # generate the x and y axis data
            x_data, y_data = self._get_plot_data(sample)

            # add lines for this sample to the plot
            fig.multi_line(x_data,
                           y_data,
                           line_width=self.line_width,
                           alpha=self.alpha,
                           color=self.colours[sample_idx],
                           legend=sample.name)
            print('Plotting finished for {}        '.format(sample.name))

        # apply axis formatting to the figure
        fig = self._axis_formatting(fig)

        if self.arguments['legend']:
            fig.legend.location = 'top_right'
            fig.legend.click_policy = 'hide'  # clicking group will hide it
            fig.legend.label_text_font = 'roboto'
            fig.legend.background_fill_alpha = 0

        if save_as:
            self.save_plot(fig, save_as)

        return fig

    def save_plot(self, fig, save_as):
        """Saves the plot to HTML.

        Args:
            fig (bokeh.plotting.Figure): A bokeh figure object.
            save_as (str): The path name to save the file as. No file extension
            should be specified.

        """
        filename = '{}.html'.format(save_as)
        output_file(filename, title=self.title)
        saved_as = bokeh_save(fig,
                              filename=filename,
                              title=self.title)
        print('Plot saved as {}'.format(saved_as))

    def _axis_formatting(self, fig):
        """Formats the axis for the line plot.

        Args:
            fig (bokeh.plotting.Figure): A bokeh Figure object.

        Returns:
            fig (bokeh.plotting.Figure): A bokeh Figure object with the axes
            formatted as required for the line plot.

        """
        # removes the x-axis ticks
        fig.xaxis.minor_tick_line_color = None
        fig.xaxis.major_tick_line_color = None
        # removes the x-axis line
        fig.xaxis.axis_line_color = None
        # limit the axis tick label bounds so they line up in the middle of each
        # base of the motif
        fig.xaxis.bounds = (0.5, len(self.motif) - 0.5)
        # force the x-axis ticks to line up with the middle of each event
        ticks = [x + 0.5 for x in range(len(self.motif))]
        fig.xaxis.ticker = FixedTicker(ticks=ticks)

        # a dictionary to map the x position to it's motif nucleotide
        label_dict = {i + 0.5: label
                      for i, label in enumerate(self.motif)}

        # creates a function that will map the x-axis label to it's motif base.
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

        for read in reads:
            # get the indices of where motif occurs within this read
            motif_idxs = read.motif_indices(self.motif)

            for i in motif_idxs:
                signal_df = read.extract_motif_signal(i)

                if signal_df.empty:
                    continue

                if self.threshold:
                    query = '{} <= {} <= {}'.format(self.threshold[0],
                                                    self.y_variable,
                                                    self.threshold[1])
                    signal_df.query(query, inplace=True)

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
