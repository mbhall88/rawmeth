import os
import pytest
from rawmeth.fast5 import Sample
from rawmeth.plotting import Plot, LinePlot


@pytest.fixture
def sample():
    """Loads in a Sample class of fast5 files contained in the test data dir.

    Returns:
        Sample: A Sample class of the fast5 files in test data dir.
    """
    return Sample('rawmeth/test/data')


@pytest.fixture
def line_plot(sample):
    """Creates the arguments for testing the line plot.

    Args:
        sample (Sample): A Sample class of the fast5 files in test data dir.

    Returns:
        (LinePlot): A figure object for the line plot

    """
    kwargs = {
        'plot_height': 750,
        'plot_width': 600,
        'title': 'Testing line plot',
        'y_axis_label': 'Signal',
        'x_axis_label': 'Motif',
        'tools': 'ypan, box_zoom, reset, save, ywheel_zoom',
        'legend': True,
        'alpha': 0.75,
        'linewidth': 3.,
        'colour_map': 'Set1',
        'threshold': (0, 2000)
    }
    y_variable = 'signal'
    motif = 'ATGGTA'
    return LinePlot(sample, motif, y_variable, **kwargs)


class TestLinePlot:

    def test_line_plot_init(self, line_plot):

        assert line_plot.alpha == 0.75
        assert line_plot.line_width == 3.
        assert line_plot.height == 750
        assert line_plot.width == 600
        assert line_plot.title == 'Testing line plot'
        assert line_plot.xlabel == 'Motif'
        assert line_plot.ylabel == 'Signal'

    def test_line_plot_setting(self, line_plot):
        line_plot.ylabel = 'ylabel'
        assert line_plot.ylabel == 'ylabel'
        line_plot.xlabel = 'xlabel'
        assert line_plot.xlabel == 'xlabel'
        line_plot.title = 'title'
        assert line_plot.title == 'title'
        line_plot.width = 550
        assert line_plot.width == 550
        line_plot.height = 980
        assert line_plot.height == 980
        line_plot.line_width = 2.
        assert line_plot.line_width == 2.
        line_plot.alpha = .25
        assert line_plot.alpha == .25
        line_plot.threshold = (0, 1500)
        assert line_plot.threshold == (0, 1500)

    def test_line_plot_save(self, line_plot):
        fname = 'rawmeth/test/data/line_plot_test.html'
        if os.path.isfile(fname):
            os.remove(fname)

        fig = line_plot.create_plot(save_as=fname.replace('.html', ''))
        assert os.path.isfile(fname)
