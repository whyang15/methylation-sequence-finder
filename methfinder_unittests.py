import os
import unittest
import tempfile
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from unittest.mock import mock_open, patch
from io import StringIO
from methFinder import find_methylase_recognition_sites
from methFinder import process_recognition_sites

class TestScriptFunctionality(unittest.TestCase):

    
    # test finding methylase recognition sites in fasta file with multiple records.
    def test_find_methylase_recognition_sites(self):
        # Set up the FASTA file path
        fasta_file = 'fasta_files/unittest.fasta'
        
        # Call the function with sample input:
        results = find_methylase_recognition_sites(fasta_file, 'TCGA', 2)

        # Verify the results:
        self.assertEqual(len(results), 2)

        # Check the first sequence:
        name, num_occurrences, positions, sequence_contexts = results[0]
        self.assertEqual(name, 'sequence1')
        self.assertEqual(num_occurrences, 2)
        self.assertEqual(positions, [1, 14])
        self.assertEqual(sequence_contexts, ['ATCGATT', 'CATCGATG'])

        # Check the second sequence
        name, num_occurrences, positions, sequence_contexts = results[1]
        self.assertEqual(name, 'sequence2')
        self.assertEqual(num_occurrences, 2)
        self.assertEqual(positions, [4, 17])
        self.assertEqual(sequence_contexts, ['GATCGATT', 'CATCGATC']) 


    # test the case where the pattern is at the beginning or end of the sequence.
    def test_find_pattern_edges(self):
        # Set up the FASTA file path
        fasta_file = 'fasta_files/unittest_edges.fasta'
    
        # Call the function with sample input:
        results = find_methylase_recognition_sites(fasta_file, 'tcga', 2)

        # Verify the results:
        self.assertEqual(len(results), 2)

        # Check the beginning sequence:
        name, num_occurrences, positions, sequence_contexts = results[0]
        # print(results)
        self.assertEqual(name, 'unittest_beginning')
        self.assertEqual(num_occurrences, 1)
        self.assertEqual(positions, [0])
        self.assertEqual(sequence_contexts, ['TCGAAG'])

        # Check the end sequence:
        name, num_occurrences, positions, sequence_contexts = results[1]
        self.assertEqual(name, 'unittest_end')
        self.assertEqual(num_occurrences, 1)
        self.assertEqual(positions, [20])
        self.assertEqual(sequence_contexts, ['CCTCGA'])


    # test the case where the script handles the case where the pattern is not found:
    def test_find_negative_ctrl(self):
        fasta_file = 'fasta_files/negative_ctrl.fasta'
        
        results = find_methylase_recognition_sites(fasta_file, 'GATC', 2)

        self.assertEqual(len(results), 1)
        
        name, num_occurrences, positions, sequence_contexts = results[0]
        self.assertEqual(name, 'negative_control_test')
        self.assertEqual(num_occurrences, 0)
        self.assertEqual(positions, [])
        self.assertEqual(sequence_contexts, [])


    # test the case where the context size is larger than the sequence length.
    def test_tiny_sequence_large_context(self):
        fasta_file = 'fasta_files/tiny_seq.fasta'
        results = find_methylase_recognition_sites(fasta_file, 'gactgact', 2)

        self.assertEqual(len(results), 1)

        name, num_occurrences, positions, sequence_contexts = results[0]
        self.assertEqual(name, 'tiny_sequence')
        self.assertEqual(num_occurrences, 0)
        self.assertEqual(positions, [])
        self.assertEqual(sequence_contexts, [])


    # test the case where the FASTA files does not exist.
    def test_nonexistent_fasta(self):
        fasta_file = 'fasta_file/positive_ctrl.fasta'

        with self.assertRaises(FileNotFoundError) as context:
            find_methylase_recognition_sites(fasta_file, 'gatc', 2)

        self.assertEqual(str(context.exception), "The input file was not found.")
    

    # test the case where the output file cannot be created.
    @patch('methFinder.open', new_callable=mock_open)
    def test_output_file_oserror(self, mock_file):

        # simulate OSError when trying to open the file
        mock_file.side_effect = OSError("Mocked OSError for testing")

        # create a temporary dummy FASTA file
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_fasta:
            temp_fasta_name = temp_fasta.name
            record = SeqRecord(Seq("ATGCATGCGATCGATCGATCGATCG"), id="Test", description="Test Sequence")
            SeqIO.write(record, temp_fasta_name, "fasta")


        # Run the function
        with self.assertRaises(SystemExit) as cm:   # Expecting a SystemExit due to 'exit(1)'
            process_recognition_sites(temp_fasta_name, 'GATC', 2, 'output.txt')

        # Clean up the temporary file:
        os.remove(temp_fasta_name)

        # Check if the exception raised was due to an error message.
        self.assertEqual(cm.exception.code, 1)

        # Ensure OSError was raised and caught
        mock_file.assert_called_once_with('output.txt', 'w')


if __name__ == '__main__':
    unittest.main()