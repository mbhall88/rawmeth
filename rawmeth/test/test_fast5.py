import pytest
from rawmeth.fast5 import Motif, Fast5


@pytest.fixture
def gatc():
    """Returns a motif class of GATC

    Returns:
        Motif: A motif of GATC

    """
    return Motif('gatc')


@pytest.fixture
def ch100():
    """Loads in a test fast5 file for testing.

    Returns:
        Fast5: A Fast5 class version of a test fast5 file.
    """
    return Fast5('rawmeth/test/data/'
                 'C4_watermang_22032017_75675_ch100_read38_strand.fast5')

@pytest.fixture
def empty_fast5():
    return Fast5('rawmeth/test/data/empty.fast5')

############################################################################
# Motif tests


class TestMotif:
    """Class to run tests for the Motif class."""

    def test_motif_init(self, gatc):
        """Tests that the motif is constructed as expected.

        Args:
            gatc (Motif): GATC motif.

        """
        assert gatc == "GATC"

    def test_motif_complement(self, gatc):
        """Tests the complement method for motif works properly.

        Args:
            gatc (Motif): GATC motif.

        """
        assert gatc.complement() == "CTAG"

    def test_motif_reverse_complement(self, gatc):
        """Tests the reverse complement method for motif works properly.

        Args:
            gatc (Motif): GATC motif.

        """
        assert gatc.reverse_complement() == "GATC"

############################################################################
# Fast5 tests


class TestFast5:
    """Class to run tests for Fast5 class."""

    def test_fast5_signal(self, ch100):
        """Tests that the signal is extracted from file correctly.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5

        """
        assert ch100.signal[2] == 508
        assert len(ch100.signal) == 40169

    def test_fast5_scale(self, ch100):
        """Tests that the nanoraw scaling is extracted from file correctly.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5

        """
        assert ch100.scale == 53.0

    def test_fast5_shift(self, ch100):
        """Tests that the nanoraw shift value is extracted from file correctly.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5

        """
        assert ch100.shift == 438.0

    def test_fast5_raw_offset(self, ch100):
        """Tests that the nanoraw read_start_rel_to_raw attribute has been
        successfully extracted.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5

        """
        assert ch100.raw_offset == 232

    def test_fast5_name(self, ch100, empty_fast5):
        """Tests that the name of the file is parsed correctly.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5
            empty_fast5 (Fast5): An empty fast5 file.

        """
        name = 'C4_watermang_22032017_75675_ch100_read38_strand.fast5'
        assert ch100.name == name
        assert empty_fast5.name == 'empty.fast5'

    def test_fast5_empty(self, ch100, empty_fast5):
        """Tests files are adjudged empty correctly.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5
            empty_fast5 (Fast5): An empty fast5 file.

        """
        assert empty_fast5.empty
        assert not ch100.empty

    def test_fast5_read_aln(self, ch100):
        """Tests read_aln from nanoraw is read in correctly.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5

        """
        assert ch100.read_aln[4655] == 'C'

    def test_fast5_genome_aln(self, ch100):
        """Tests genome_aln from nanoraw is read in correctly.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5

        """
        assert ch100.genome_aln[468] == 'T'

    def test_fast5_events(self, ch100):
        """Tests that the nanoraw events table is extracted correctly.

        Args:
            ch100 (Fast5): A Fast5 class structure form of the test file
            C4_watermang_22032017_75675_ch100_read38_strand.fast5

        """
        a, b, c, d, e = ch100.events[4399]
        assert a == 1.4177897574123988
        assert b == 0.3021273589046265
        assert c == 38059
        assert d == 7
        assert e == 'C'

    def test_fast5_motif_indices(self, ch100):
        """Tests that the motif index extraction is working properly.

        Args:
            ch100 ():
        """
        assert ch100.motif_indices('GATC') == [(692, 696), (701, 705),
                                               (749, 753), (866, 870),
                                               (878, 882), (954, 958),
                                               (1505, 1509), (1583, 1587),
                                               (1671, 1675), (1730, 1734),
                                               (1833, 1837), (2020, 2024),
                                               (2403, 2407), (2466, 2470),
                                               (2613, 2617), (2781, 2785),
                                               (3390, 3394), (3650, 3654),
                                               (3696, 3700), (3764, 3768),
                                               (4002, 4006), (4014, 4018),
                                               (4386, 4390)]
        assert ch100.motif_indices('TCAACGAAATC') == [(14, 25)]

    # TODO: finish test cases