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

    def __init__(self, path, limit=None):

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
        # type: (str) -> Fast5
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
        # type: (str) -> None

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
        # type: () -> bool
        """Function to determine if a fast5 file has been basecalled.
        This is based on whether it has the 'Analyses' group.

        :return: True or False
        """
        return 'Analyses' in self.f.keys()

    def _is_corrected(self):
        # type: () -> bool
        """Function to check whether the fast5 file has been resquiggled by
        nanoraw.

        :return: Whether or not there is corrected data.
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
        # type: (str, object) -> None
        """Function to walk the fast5 internal directory and extract the
        necessary information.

        :param name: The name of the current group
        :param obj: The information attached to the group
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
        # type: (Motif) -> list[tuple[int, int]]
        """Will find motif occurrances (overlapping) in ungapped sequence

        :param motif: The DNA motif to find indices for. i.e 'GATC'
        :return: A list of tuples containing the start and end index for that
        motif in the sequence (events['base']).
        """
        if not isinstance(motif, Motif):
            motif = Motif(motif)

        seq = ''.join(self.events['base'])
        return [m.span(1)
                for m in re.finditer(r'(?=({}))'.format(motif.regex()), seq)]

    def extract_motif_signal(self, idx):
        # type: (tuple[int, int]) -> pd.DataFrame
        """For a given start/end index for a motif, will extract the raw signal
        associated with each base in the motif.

        :param idx: tuple containing the start and end index within in the raw
        signal array that the motif maps to. (start, end)
        :return: A Pandas DataFrame. Each row has the raw signal, the base that
        signal matches to, and the median normalised raw signal.
        """
        starts = []
        lengths = []
        bases = []
        motif_events = self.events[slice(*idx)]
        for e in motif_events:
            s, l, b = list(e)[2:]
            starts.append(s)
            lengths.append(l)
            bases.append(b)

        d = self._extract_raw_signal(starts, lengths, bases)
        df = pd.DataFrame.from_dict(d)
        df['motif'] = ''.join(bases)
        return df

    def get_motif_signals(self, motif):
        # type: (Motif) -> pd.DataFrame
        """Will return the raw signals associated with all occurrences of a
         given motif.

        :param motif: DNA motif of interest. i.e 'GATC'
        :return: A Pandas DataFrame. Each row has the raw signal, the base that
        signal matches to, and the median normalised raw signal.
        """
        if not isinstance(motif, Motif):
            motif = Motif(motif)

        idxs = self.motif_indices(motif)
        all_dfs = [df
                   for df in map(self.extract_motif_signal, idxs)
                   if not df.empty]

        return pd.concat(all_dfs) if all_dfs else pd.DataFrame()

    def _extract_raw_signal(self, starts, lengths, bases):
        # type: (list[int], list[int], list[str]) -> dict
        """Maps the information for each event onto the raw signal array and
        extracts it.

        :param starts: List of the indices (for a motif) that denote the index
         for the raw signal at the beginning of an event.
        :param lengths: List of lengths of each event in motif.
        :param bases: Bases that make up the motif of interest.
        :return: A Pandas DataFrame. Each row has the raw signal, the base that
        signal matches to, and the median normalised raw signal.
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

            # assign a base to each raw signal and create a dataframe
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
        idxs = self.motif_indices(motif)
        rows_list = []
        for i in idxs:
            motif_events = self.events[slice(*i)]
            for pos, row in enumerate(motif_events):
                l, b = list(row)[3:]
                rows_list.append({
                    'length': l,
                    'base': b,
                    'pos': pos,
                    'motif': motif
                })
        return pd.DataFrame(rows_list)


def flatten_list(xs):
    # type: (list) -> list
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
