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
from itertools import chain
import h5py as h5
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
        d = {
            'A': 'T', 'T': 'A',
            'C': 'G', 'G': 'C',
            'W': 'S', 'S': 'W',
            'M': 'K', 'K': 'M',
            'R': 'Y', 'Y': 'R',
            'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D',
            'N': 'N'
        }
        return Motif("".join([d[c] for c in self]))

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
        d = {
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
        return ''.join([d[c] for c in self])


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

    @staticmethod
    def _load_f5(filename):
        """Loads a fast5 file in and initiates it as a Fast5 class.

        Args:
            filename (str): The path to the fast5 file to load in.

        Returns:
            f5 (Fast5): Returns the Fast5 representation of the file, or, if
            the file hasn't been basecalled or nanoraw corrected, will return
            None.

        """
        f5 = Fast5(filename)
        if not f5.empty:
            return f5

    def get_motif_signal(self, motifs):
        """Constructs a dataframe of all raw signals for a motif within sample.

        Extracts the raw signal corresponding to the given motif(s) for all
        fast5 files in the sample and returns them in a single dataframe.

        Args:
            motifs (str | list[str]): A DNA motif or a list of DNA motifs.

        Returns:
            df (pd.DataFrame): A dataframe with all raw signals for given
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

        df = pd.concat(motif_dfs)
        df['sample'] = self.basename
        return df


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
            with h5.File(path, 'r') as f:
                self.f = f
                if self._is_corrected():
                    self.f.visititems(self.extract_fast5_info)
                    self.empty = False

        except IOError as e:
            msg = "Unable to open file (File signature not found)"
            if msg == e.message:
                self.empty = True
            else:
                raise e

    def _is_basecalled(self):
        """Checks if a fast5 file has been basecalled.

        This is based on whether it has the 'Analyses' group.

        Returns:
            bool

        """
        return 'Analyses' in self.f.keys()

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
        for k in self.f['Analyses'].keys():
            if nanoraw_re.match(k):
                if not self.f['Analyses/' + k].values():
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
            self.signal = self.f[name].value
        elif events_re.match(name):
            self.events = self.f[name].value
            self.raw_offset = self.f[name].attrs['read_start_rel_to_raw']
        elif read_segs_re.match(name):
            self.read_segs = self.f[name].value
        elif read_aln_re.match(name):
            self.read_aln = self.f[name].value
        elif genome_aln_re.match(name):
            self.genome_aln = self.f[name].value
        elif attrs_re.match(name):
            self.scale = self.f[name].attrs['scale']
            self.shift = self.f[name].attrs['shift']

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
        for e in motif_events:
            s, length, b = list(e)[2:]
            starts.append(s)
            lengths.append(length)
            bases.append(b)

        d = self._extract_raw_signal(starts, lengths, bases)
        df = pd.DataFrame.from_dict(d)
        df['motif'] = ''.join(bases)
        return df

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
                length, b = list(row)[3:]
                rows_list.append({
                    'length': length,
                    'base': b,
                    'pos': pos,
                    'motif': motif
                })
        return pd.DataFrame(rows_list)


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
