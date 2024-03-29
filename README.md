# data: database, and test sequence

The `dog_breeds.fa` file contains the sequence in our database, to compare against.
The `mystery.fa` sequence contains a sequence from an unknown dog type, and we want to identify what is the closest breed to it in our database.
The `shamima/main.py` file contains Python code for comparison.
The `shamima/results.txt` file contains the most compatible breed and alignment that is nearly same, with similarity score.

# Phylogeny

If you want to make a phylogeny, the simplest way to do it is to feed the file to an online webservice, such as https://www.ebi.ac.uk/Tools/phylogeny/simple_phylogeny/
export a file, then read it in using the Phylo submodule.

You can also compute a phylogenetic tree via the following two steps:
- construct a distance matrix
- construct a tree using a construction algorithm
See the following page for examples: https://biopython.org/docs/1.75/api/Bio.Phylo.TreeConstruction.html

# Code Overview

Function: read_dna_files

Description:
This function reads a DNA file and extracts DNA sequences along with their corresponding dog breeds.

Parameters:
- filepath: str
  The path to the DNA file containing dog breed information.

Returns:
- If only one breed is present in the file, the DNA sequence for that breed is returned as a string.
- If multiple breeds are present, a list of lists is returned, where each inner list contains a dog breed and its corresponding DNA sequence.

Function: compare_alignments

Description:
This function performs a global pairwise alignment between a known breed DNA sequence and a mystery breed DNA sequence.

Parameters:
- mystery_breed_dna_sequence: str
  The DNA sequence of the mystery breed.
- known_breed_dna_sequence: str
  The DNA sequence of the known breed for comparison.

Returns:
- alignments: list
  A list of alignments representing the pairwise alignments between the two DNA sequences.

Function: calculate

Description:
This function calculates the compatibility between a mystery DNA sequence and known dog breeds' DNA sequences.

Parameters:
- dog_breed_with_dna: list
  A list of lists containing dog breeds and their corresponding DNA sequences.
- mystery_dna_sequence: str
  The DNA sequence of the mystery breed.

Returns:
- None

Output:
- The function writes the results of the compatibility calculation to a file named "results.txt" in the current directory.
- Additionally, it prints the most compatible breed, the maximum alignment score, and the alignment details for the most compatible sequence.

Main Execution:
- Reads known dog breeds DNA file and mystery dog breed DNA file.
- Calls the calculate function to perform compatibility calculations.

Test Execution:
- Used a known dog breed as a mystery dog breed (mystery_test.fa).
- Run the test_main.py file.
