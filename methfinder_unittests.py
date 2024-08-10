import unittest
from io import StringIO
import sys
from methFinder import main

class TestScriptFunctionality(unittest.TestCase):
    def setUp(self):
        # Save the original stdout to restore it later
        sys.original_stdout = sys.stdout

    def tearDown(self):
        # Restore the original stdout
        sys.stdout = self.original_stdout

    def test_pattern_found(self):
        # Capture the output
        captured_output = StringIO()
        sys.stdout = captured_output

        # Run the script with a FASTA file containing the pattern
        main(['--input', 'fasta_files/dam_test.fasta', '--find', 'GATC', '--numbases', '2'])

        # Assert the output contains the expected sequence contexts
        self.assertIn('GCGATCAA', captured_output.getvalue())
    