"""This is a module to provide functions to handle basic data structures for
fast5 files and collections of fast5 files for a sample/experiment. There is
also a class to handle DNA sequence motifs (allows for use of ambiguous bases).

This module is designed with the idea of using the dataframes produced to
plot and explore the raw signal associated with given DNA motifs.

"""

from __future__ import print_function
import re
import glob
import os
import sys
from itertools import chain
from collections import Counter
import h5py as h5
import pandas as pd


LENGTH_FILTER = 100


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

    def __init__(self, path, limit=None, name=None):
        """Initiates the Sample class from a directory.

        Args:
            path (str): Path to the directory of fast5 files.
            limit (int): Used if only wanting to load in a certain number of
             files. i.e if you have 100,000 files but just want to explore a
             smaller subset.
            name (str): Name to designate for the sample. Will default to the
            name of the directory the files are contained within.

        """
        if os.path.isdir(path):
            self.path = os.path.normpath(path)  # remove any trailing '/'
            self.basename = name or os.path.basename(self.path)
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
        self.files = []
        num_file_paths = len(self.file_paths)

        for i, filename in enumerate(self.file_paths):

            fast5 = self._load_f5(filename)

            if fast5:
                self.files.append(fast5)

            perc_complete = float(i) / num_file_paths * 100
            print('{}% of files loaded for {}'
                  '        '.format(perc_complete, self.basename), end='\r')
            sys.stdout.flush()

        print('All files loaded for {}            '.format(self.basename))

    def __getitem__(self, item):
        """Allows for indexing the Sample object.

        Args:
            item (int): Index

        Returns:
            (Fast5): The Fast5 object within self.files list corresponding to
            the given index.

        """
        return self.files[item]

    def __iter__(self):
        """When iterating on Sample object, will iterate over fast5 files."""
        return iter(self.files)

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

    @property
    def size(self):
        """Returns the number of files in the Sample."""
        return len(self.files)

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

    def motif_counts(self, motif, pretty_print=False):
        """Returns the counts for all possible variations of a given motif.

        Examples:
            >>> import re
            >>> from collections import Counter
            >>> seq = 'GATTGAGTGAGTGAGTGAATGAAT'
            >>> my_re = re.compile(r'GA[ACGT]T')
            >>> matches = my_re.findall(seq)
            >>> Counter(matches)
            Counter({'GAGT': 3, 'GAAT': 2, 'GATT': 1})

        Args:
            motif (Motif): A DNA motif to count. Can include ambiguous bases.
            pretty_print (bool): Whether to print the counts nicely before
            returning the counter.

        Returns:
            (Counter): A subclass of a dictionary. The keys are the motif and
            the values are the counts for that motif.

        """
        # get list of all counters
        counts = [fast5.motif_counts(motif) for fast5 in self]

        # sum the list of counters by combining them
        combined_counts = sum(counts, Counter())

        if pretty_print:
            pretty_print_counts(combined_counts)

        return combined_counts


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

        # fixme: fix this disgusting mess!!!
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
        motif = Motif(motif)

        seq = ''.join(self.events['base'])
        return [m.span(1)
                for m in re.finditer(r'(?=({}))'.format(motif.regex()), seq)]

    def motif_counts(self, motif, pretty_print=False):
        """Returns the counts for all possible variations of a given motif.

        Examples:
            >>> import re
            >>> from collections import Counter
            >>> seq = 'GATTGAGTGAGTGAGTGAATGAAT'
            >>> my_re = re.compile(r'GA[ACGT]T')
            >>> matches = my_re.findall(seq)
            >>> Counter(matches)
            Counter({'GAGT': 3, 'GAAT': 2, 'GATT': 1})

        Args:
            motif (Motif): A DNA motif to count. Can include ambiguous bases.
            pretty_print (bool): Whether to print the counts nicely before
            returning the counter.

        Returns:
            (Counter): A subclass of a dictionary. The keys are the motif and
            the values are the counts for that motif.

        """
        motif = Motif(motif)

        seq = ''.join(self.events['base'])
        matches = re.findall(r'(?=({}))'.format(motif.regex()), seq)

        counts = Counter(matches)

        if pretty_print:
            pretty_print_counts(counts)

        return counts

    def extract_motif_signal(self, idx):
        """For a given start/end index for a motif, will extract the raw signal
        associated with each base in the motif.

        Args:
            idx (tuple[int, int]): tuple containing the start and end index
            within in the raw signal array that the motif maps to. (start, end).

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

            # todo: make this a little more robust
            # filter out motifs with events longer than 100
            if length > LENGTH_FILTER:
                return pd.DataFrame()

            starts.append(start)
            lengths.append(length)
            bases.append(base)

        signal_dict = self._extract_raw_signal(starts, lengths, bases)

        # if the dictionay is empty, return an empty dataframe
        if not signal_dict:
            return pd.DataFrame()

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


def flatten_list(xs):
    """Completely flattens a list to give a single list.

    Args:
        xs (list): A nested list of any level of nesting.

    Returns:
        list: A list with no nesting.

    """
    return list(chain.from_iterable(xs))


def _is_list_empty(xs):
    """Determines if a list is truly empty.

    Goes through a given list recursively and checks whether all sequence
    elements inside are empty.

    Example:
        xs = [{}, [], set(), '']
        _is_list_empty(xs)  # True

    Args:
        xs (list): List of anything.

    Returns:
        bool: Whether or not list is truly empty.

    """
    try:
        return all(_is_list_empty(x) for x in xs)
    except TypeError:
        return False


def pretty_print_counts(counts):
    """Prints the given counts nicely, with the most common on top and the
    least common on the bottom.

    Args:
        counts (Counter): An instance of a Counter, which is a dictionary of
        counts for a motif.

    """
    sorted_counts = counts.most_common()
    motif_len = len(sorted_counts[0][0])
    motif_header = 'Motif'

    # add some 'padding' (spaces) to header if motif is long
    if len(motif_header) < motif_len:
        motif_header += ' ' * (motif_len - len(motif_header))

    # print header and separator
    header = '{}    Counts'.format(motif_header)
    print(header)
    print('-' * len(header))

    # loop through and print motif and counts
    for motif, count in sorted_counts:
        print('{}\t{}'.format(motif, count))
