"""This is a module to provide functions to handle basic data structures for
fast5 files and collections of fast5 files for a sample/experiment. There is
also a class to handle DNA sequence motifs (allows for use of ambiguous bases).

This module is designed with the idea of using the dataframes produced to
plot and explore the raw signal associated with given DNA motifs.

"""
# TODO: methods combining Sample dataframes

from __future__ import print_function
import re
import glob
import os
from itertools import chain, groupby
import h5py as h5
import matplotlib.pyplot as plt
import pandas as pd


class Motif(str):
    """Class that is used to store a DNA motif.

    Acts effectively like the str class but with some additional methods for
    DNA manipulation and ambiguous bases.

    Ambiguous bases are based on the IUPAC naming convention.

    """

    def __new__(cls, s):
        """Creates a new DNA motif. Converts all characters to uppercase.

        Args:
            s (str): A DNA motif.

        Returns:
            str: Uppercase version of the given str.

        """
        return super(Motif, cls).__new__(cls, s.upper())

    def complement(self):
        """Returns the DNA complement of the motif."""
        dna_table = {
            'A': 'T', 'T': 'A',
            'C': 'G', 'G': 'C',
            'W': 'S', 'S': 'W',
            'M': 'K', 'K': 'M',
            'R': 'Y', 'Y': 'R',
            'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D',
            'N': 'N'
        }
        return Motif("".join([dna_table[c] for c in self]))

    def reverse_complement(self):
        """Returns the DNA reverse complement of the motif."""
        return Motif(self.complement()[::-1])

    def regex(self):
        """Constructs a skeleton that can be used to create a regular
        expression for the motif.

        Method doesn't return a regular expression class, but this can simply
        be done with the following:

        import re
        motif = Motif('GWAD')
        my_regex = re.compile(r'{}'.format(motif.regex()))

        Returns:
            str: A string representation for the motif as a regular expression.

        """
        dna_table = {
            'A': 'A',
            'C': 'C',
            'G': 'G',
            'T': 'T',
            'W': '[AT]',
            'S': '[CG]',
            'M': '[AC]',
            'K': '[GT]',
            'R': '[AG]',
            'Y': '[CT]',
            'B': '[^A]',
            'D': '[^C]',
            'H': '[^G]',
            'V': '[^T]',
            'N': '.'
        }
        return ''.join([dna_table[c] for c in self])


class Sample(object):
    """This is a class that holds all of the Fast5 objects in a single class.
    You provide a path to a directory containing fast5 files belonging to a
    sample/experiment and can then extract information for all of the reads in
    that sample.

    """

    def __init__(self, path, limit=None):
        """Initiates the Sample class from a directory.

        Args:
            path (str): Path to the directory of fast5 files.
            limit (int): Used if only wanting to load in a certain number of
             files. i.e if you have 100,000 files but just want to explore a
             smaller subset.

        """
        if os.path.isdir(path):
            self.path = os.path.normpath(path)  # remove any trailing '/'
            self.basename = os.path.basename(self.path)
        else:
            raise Exception(
                "{} does not exist. Provide a valid directory.".format(path)
            )

        self.file_paths = glob.glob(os.path.join(self.path, '*.fast5'))[:limit]

        if not self.file_paths:
            raise Exception(
                "No '.fast5' files could be found in {}.".format(self.path)
            )

        # load fast5 files into Fast5 objects
        self.files = filter(None, map(self._load_f5, self.file_paths))

    def __iter__(self):
        """When iterating on Sample object, will iterate over fast5 files."""
        for fast5 in self.files:
            yield fast5

    @staticmethod
    def _load_f5(filename):
        """Loads a fast5 file in and initiates it as a Fast5 class.

        Args:
            filename (str): The path to the fast5 file to load in.

        Returns:
            fast5 (Fast5): Returns the Fast5 representation of the file, or, if
            the file hasn't been basecalled or nanoraw corrected, will return
            None.

        """
        fast5 = Fast5(filename)
        if not fast5.empty:
            return fast5

    @property
    def name(self):
        """Returns the current value for the name of the sample."""
        return self.basename

    @name.setter
    def name(self, name):
        """Sets the name of the sample.

        Args:
            name (str): Name to assign to sample.
        """
        self.basename = name

    def get_motif_signal(self, motifs):
        """Constructs a dataframe of all raw signals for a motif within sample.

        Extracts the raw signal corresponding to the given motif(s) for all
        fast5 files in the sample and returns them in a single dataframe.

        Args:
            motifs (str | list[str]): A DNA motif or a list of DNA motifs.

        Returns:
            master_df (pd.DataFrame): A dataframe with all raw signals for given
            motif(s) in sample.

        """
        # make sure motif is of class Motif and is in list form.
        if isinstance(motifs, str):
            motifs = [Motif(motifs)]
        else:
            motifs = [Motif(m) for m in motifs]

        motif_dfs = []
        for motif in motifs:
            dfs = [f5.get_motif_signals(motif) for f5 in self.files]
            motif_df = pd.concat(dfs)
            motif_dfs.append(motif_df)

        master_df = pd.concat(motif_dfs)
        master_df['sample'] = self.name
        return master_df

    def line_plot(self, motif, against=None, yaxis='signal', alpha=None,
                  linewidth=None, colour_map='Set1', save=None):
        """Produces a line plot of the raw signal events related to the given
         motif.

        Args:
            motif (Motif): The motif to plot signal for.
            against (list[Sample]): A list of Samples to compare against.
            yaxis (str): Variable to plot on the y-axis.
            alpha (float): Transparency of the lines. Must be in the range
            [0, 1]
            linewidth (float): Width of the lines.
            colour_map (str): Matplotlib colour map to use. For list of names -
            https://matplotlib.org/examples/color/colormaps_reference.html
            save (str): Filename to save plot as. If left as None, file will
            not be saved.

        """
        colours = plt.get_cmap(colour_map).colors

        # make sure anything wanting to plot against is in a list form.
        if not isinstance(against, list):
            if against:
                against = [against]
            else:
                against = []

        against.append(self)
        # loop through sample files and plot a line for each of their occurrence
        # of the motif.
        for idx, sample in enumerate(against):
            for fast5 in sample.files:
                motif_idxs = fast5.motif_indices(motif)
                for i in motif_idxs:
                    signal_df = fast5.extract_motif_signal(i)
                    x = generate_line_plot_xs(signal_df['pos'])
                    y = signal_df[yaxis]
                    plt.plot(x, y, linewidth=linewidth, alpha=alpha,
                             color=colours[idx])

        if save:
            plt.savefig(save, dpi=300)

        plt.show()
        plt.close()



class Fast5(object):
    """
    Holds the information extracted from a Fast5 sequence file. Also provides
    methods for extracting information about motifs and their associated raw
    signal from the file.

    """

    def __init__(self, path):
        """Initiates the Fast5 class from a file path.

        Args:
            path (str): Path to the fast5 to load in.

        """
        self.read_aln = None
        self.genome_aln = None
        self.signal = None
        self.events = None
        self.raw_offset = None
        self.read_segs = None
        self.scale = None
        self.shift = None
        self.name = path.split('/')[-1]
        self.empty = True

        # open the fast5 file and extract the required information, then close.
        try:
            with h5.File(path, 'r') as file_:
                self.file_ = file_
                if self._is_corrected():
                    self.file_.visititems(self.extract_fast5_info)
                    self.empty = False

        except IOError as err:
            msg = "Unable to open file (File signature not found)"
            if msg == err.message:
                self.empty = True
            else:
                raise err

    def _is_basecalled(self):
        """Checks if a fast5 file has been basecalled.

        This is based on whether it has the 'Analyses' group.

        Returns:
            bool

        """
        return 'Analyses' in self.file_.keys()

    def _is_corrected(self):
        """Checks whether the fast5 file has been resquiggled by nanoraw.

        Checks that there is actually information nested under the
        'RawGenomeCorrected' group.

        Returns:
            bool: Whether or not there is corrected data.

        """
        if not self._is_basecalled():
            return False

        nanoraw_re = re.compile(r'RawGenomeCorrected*')
        for k in self.file_['Analyses'].keys():
            if nanoraw_re.match(k):
                if not self.file_['Analyses/' + k].values():
                    return False
                return True

    def extract_fast5_info(self, name, obj):
        """Function to walk the fast5 internal directory and extract the
        necessary information.

        Args:
            name (str): The name of the current group
            obj (object): The information attached to the group

        """
        # create regular expression for searching for data
        skel = r'Analyses/RawGenomeCorrected_\d+/BaseCalled_template'
        attrs_re = re.compile(r'{}$'.format(skel))
        events_re = re.compile(r'{}/Events$'.format(skel))
        read_segs_re = re.compile(r'{}/Alignment/read_segments$'.format(skel))
        read_aln_re = re.compile(r'{}/Alignment/read_alignment$'.format(skel))
        genome_aln_re = re.compile(
            r'{}/Alignment/genome_alignment$'.format(skel))
        signal_re = re.compile(r'Raw/Reads/Read_\d+/Signal$')

        # TODO: fix this disgusting mess!!!
        if signal_re.match(name):
            self.signal = self.file_[name].value
        elif events_re.match(name):
            self.events = self.file_[name].value
            self.raw_offset = self.file_[name].attrs['read_start_rel_to_raw']
        elif read_segs_re.match(name):
            self.read_segs = self.file_[name].value
        elif read_aln_re.match(name):
            self.read_aln = self.file_[name].value
        elif genome_aln_re.match(name):
            self.genome_aln = self.file_[name].value
        elif attrs_re.match(name):
            self.scale = self.file_[name].attrs['scale']
            self.shift = self.file_[name].attrs['shift']

    def motif_indices(self, motif):
        """Will find motif occurrances (overlapping) in ungapped sequence

        Args:
        motif (Motif): The DNA motif to find indices for. i.e 'GATC'

        Returns:
            list[tuple[int, int]]: A list of tuples containing the start and
            end index for that motif in the sequence (events['base']).

        """
        if not isinstance(motif, Motif):
            motif = Motif(motif)

        seq = ''.join(self.events['base'])
        return [m.span(1)
                for m in re.finditer(r'(?=({}))'.format(motif.regex()), seq)]

    def extract_motif_signal(self, idx):
        """For a given start/end index for a motif, will extract the raw signal
        associated with each base in the motif.

        Args:
            idx (tuple[int, int]): tuple containing the start and end index
            within in the raw signal array that the motif maps to. (start, end)

        Returns:
            pd.DataFrame: Each row has the raw signal, the base that
            signal matches to, the median normalised raw signal, and the
            position within the motif for that event.

        """
        starts = []
        lengths = []
        bases = []
        motif_events = self.events[slice(*idx)]
        for event in motif_events:
            start, length, base = list(event)[2:]
            starts.append(start)
            lengths.append(length)
            bases.append(base)

        signal_dict = self._extract_raw_signal(starts, lengths, bases)
        signal_df = pd.DataFrame.from_dict(signal_dict)
        signal_df['motif'] = ''.join(bases)
        return signal_df

    def get_motif_signals(self, motif):
        """Will return the raw signals associated with all occurrences of a
         given motif.

        Args:
        motif (Motif): DNA motif of interest. i.e 'GATC'.

        Returns:
            pd.DataFrame: Each row has the raw signal, the base that
            signal matches to, the median normalised raw signal, and the
            position within the motif for that event.

        """
        if not isinstance(motif, Motif):
            motif = Motif(motif)

        idxs = self.motif_indices(motif)
        all_dfs = [df
                   for df in map(self.extract_motif_signal, idxs)
                   if not df.empty]

        return pd.concat(all_dfs) if all_dfs else pd.DataFrame()

    def _extract_raw_signal(self, starts, lengths, bases):
        """Maps the information for each event onto the raw signal array and
        extracts it.

        Args:
            starts (list[int]): List of the indices (for a motif) that denote
            the index for the raw signal at the beginning of an event.
            lengths (list[int]): List of lengths of each event in motif.
            bases (list[str]): Bases that make up the motif of interest.

        Returns:
            dict: A dictionary. Each 'row' has the raw signal, the base that
        signal matches to, the median normalised raw signal, and the position
        within the motif for that event.

        """
        try:
            # create list of motif positions as strings which will be used
            # for investigating positions within a motif. If the same base
            # occurs multiple times in a motif, just investigating based on
            # base label will grab data for multiple positions within a motif
            positions = map(str, range(len(bases)))

            # make sure the starts index doesn't exceed the signals index
            idx_of_last_signal = starts[-1] + self.raw_offset + lengths[-1]
            assert idx_of_last_signal < len(self.signal), \
                "Index for starts is outside the domain of signal"

            # get the raw signals that map to the motif
            start_sig_idx = starts[0] + self.raw_offset
            end_sig_idx = starts[-1] + lengths[-1] + self.raw_offset
            sigs = self.signal[start_sig_idx:end_sig_idx]
            labels = flatten_list([list(b) * l
                                   for l, b in zip(lengths, bases)])
            pos = flatten_list([[p] * l
                                for l, p in zip(lengths, positions)])

            # assign a base to each raw signal and create a dictionary
            return {
                'signal': sigs,
                'base': labels,
                'signal_norm': [(x - self.shift) / self.scale for x in sigs],
                'pos': pos
            }

        # catch reads where the start index exceeds the length of the signal
        except AssertionError:
            # reads_skipped.append(f)
            return {}

    def get_motif_lengths(self, motif):
        """Get a dataframe of all the event lengths for each base in a motif
        across the read.

        Args:
            motif (Motif): A motif whose base event lengths to extract.

        Returns:
            pd.DataFrame: A dataframe with each row containing a length (int),
            base (str), position (str), and motif (str).

        """
        idxs = self.motif_indices(motif)
        rows_list = []
        for i in idxs:
            motif_events = self.events[slice(*i)]
            for pos, row in enumerate(motif_events):
                length, base = list(row)[3:]
                rows_list.append({
                    'length': length,
                    'base': base,
                    'pos': pos,
                    'motif': motif
                })
        return pd.DataFrame(rows_list)

    def line_plot(self, motif, against=None, yaxis='signal', alpha=None,
                  linewidth=None, colour_map='Set1', save=None):
        """Produces a line plot of the raw signal events related to the given
         motif.

        Args:
            motif (Motif): The motif to plot signal for.
            against (list[Fast5]): A list of Fast5 files to compare against.
            yaxis (str): Variable to plot on the y-axis.
            alpha (float): Transparency of the lines. Must be in the range
            [0, 1]
            linewidth (float): Width of the lines.
            colour_map (str): Matplotlib colour map to use. For list of names -
            https://matplotlib.org/examples/color/colormaps_reference.html
            save (str): Filename to save plot as. If left as None, file will
            not be saved.

        """
        colours = plt.get_cmap(colour_map).colors

        # make sure anything wanting to plot against is in a list form.
        if not isinstance(against, list):
            if against:
                against = [against]
            else:
                against = []

        against.append(self)
        # loop through fast5 files and plot a line for each of their occurrence
        # of the motif.
        for idx, fast5 in enumerate(against):
            motif_idxs = fast5.motif_indices(motif)
            for i in motif_idxs:
                signal_df = fast5.extract_motif_signal(i)
                x = generate_line_plot_xs(signal_df['pos'])
                y = signal_df[yaxis]
                plt.plot(x, y, linewidth=linewidth, alpha=alpha,
                         color=colours[idx])

        if save:
            plt.savefig(save, dpi=300)

        plt.show()
        plt.close()

def flatten_list(xs):
    """Completely flattens a list to give a single list.

    Args:
        xs (list): A nested list of any level of nesting.

    Returns:
        list: A list with no nesting.

    """
    return list(chain.from_iterable(xs))


def is_list_empty(xs):
    """Determines if a list is truly empty.

    Goes through a given list recursively and checks whether all sequence
    elements inside are empty.

    Example:
        xs = [{}, [], set(), '']
        is_list_empty(xs)  # True

    Args:
        xs (list): List of anything.

    Returns:
        bool: Whether or not list is truly empty.

    """
    try:
        return all(is_list_empty(x) for x in xs)
    except TypeError:
        return False


def generate_line_plot_xs(xs):
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
        generate_line_plot_xs(xs)
        >>> [0.0, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 0.0, 0.5]

    Args:
        xs (list): A list of integers or floats.

    Returns:
        x_coords (list[float]): A list of the x-coordinates to use to
        'squash` the given events onto the same x-axis domain.

    """
    x_coords = []
    xs_generator = groupby(xs)
    for key, gen in xs_generator:
        event_len = sum(1 for _ in gen)
        pos_as_perc = [(float(pos) / event_len) + int(key)
                       for pos in range(event_len)]
        x_coords.extend(pos_as_perc)
    return x_coords
