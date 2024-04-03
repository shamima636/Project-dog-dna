# Import necessary modules from Biopython
from Bio import SeqIO  # Used to parse sequences from files
from Bio.Phylo.TreeConstruction import DistanceCalculator, \
    DistanceTreeConstructor  # Used to calculate distances between sequences and construct trees
from Bio.Phylo import draw  # Used to draw phylogenetic trees
from Bio.SeqRecord import SeqRecord  # Used to represent individual sequences as SeqRecord objects
from Bio.Align import MultipleSeqAlignment  # Used to create multiple sequence alignments


# Define a function to read DNA sequences from a file
def read_dna_files(filepath):
    # Initialize an empty dictionary to store breed sequences
    breeds = {}
    # Iterate over each sequence record in the file
    for record in SeqIO.parse(filepath, "fasta"):
        # Extract breed information from the description of the sequence
        breed = record.description.split("[")[7].split("=")[1][:-2]
        # Extract the DNA sequence as a string
        sequence = str(record.seq)
        # Store the breed and its corresponding DNA sequence in the dictionary
        breeds[breed] = sequence
    # Return the dictionary containing breed sequences
    return breeds


# Define a function to calculate the similarity between a mystery DNA sequence and known dog breeds
def calculate(dog_breeds, mystery_dna_sequence):
    # Initialize variables to store the best match information
    matching_breed = None
    max_score = 0
    closest_sequence = None
    closest_difference = None
    closest_alignment = None

    # Initialize a list to store distances between the mystery sequence and known breeds
    distances = []

    # Iterate over each known dog breed and its DNA sequence
    for breed, dna_sequence in dog_breeds.items():
        # Create SeqRecord objects for the known breed sequence and the mystery DNA sequence
        record1 = SeqRecord(dna_sequence, id="Seq1")
        record2 = SeqRecord(mystery_dna_sequence, id="Seq2")

        # Create a MultipleSeqAlignment object containing the two sequences
        alignment = MultipleSeqAlignment([record1, record2])

        # Create a DistanceCalculator object with 'identity' as the method
        calculator = DistanceCalculator('identity')

        # Calculate the distance matrix between the sequences in the alignment
        matrix = calculator.get_distance(alignment)

        # Compute the similarity score (1 - identity distance) between the sequences
        score = 1 - matrix[0][1]

        # Append the score to the distances list
        distances.append((breed, score, matrix))

        # Update the best match if the current score is higher than the maximum score
        if score > max_score:
            max_score = score
            matching_breed = breed
            closest_sequence = dna_sequence
            closest_alignment = alignment
            closest_difference = sum(1 for a, b in zip(mystery_dna_sequence, dna_sequence) if a != b)

    # Write the distances to a file named "distances.txt"
    with open("results.txt", "w") as f:
        for dist in distances:
            f.write(f"{dist[0]}: {dist[1]}\n")

    # Return the closest DNA sequence and its difference
    return closest_sequence, closest_difference, closest_alignment, matrix, matching_breed


# Define a function to generate a phylogenetic tree
def generate_tree(dog_breeds):
    # Initialize a list to store SeqRecord objects for each breed
    seq_records = []

    # Iterate over each known dog breed and its DNA sequence
    for breed, dna_sequence in dog_breeds.items():
        # Create a SeqRecord object for the breed sequence
        record = SeqRecord(dna_sequence, id=breed)
        # Append the SeqRecord to the list
        seq_records.append(record)

    # Create a MultipleSeqAlignment object containing all sequences
    alignment = MultipleSeqAlignment(seq_records)

    # Create a DistanceCalculator object with 'identity' as the method
    calculator = DistanceCalculator('identity')

    # Calculate the distance matrix between the sequences in the alignment
    matrix = calculator.get_distance(alignment)

    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()

    # Build the tree using the distance matrix
    tree = constructor.upgma(matrix)

    # Draw the phylogenetic tree"
    draw(tree, branch_labels=lambda c: round(c.branch_length * 100, 2))


# Define a function to run the comparison between mystery and known DNA sequences
def run(test_file):
    # Read the known dog breeds DNA sequences from the provided file
    dog_breeds = read_dna_files("../dog_breeds.fa")

    # Read the mystery dog breed DNA sequence from the test file
    mystery_dna_sequence = list(read_dna_files(test_file).values())[0]

    # Calculate the closest matching breed for the mystery DNA sequence
    closest_sequence, closest_difference, closest_alignment, distance_matrix, matching_breed = calculate(dog_breeds, mystery_dna_sequence)

    # Generate and save the phylogenetic tree
    generate_tree(dog_breeds)

    # Return the closest DNA sequence and its difference
    return closest_sequence, closest_difference, closest_alignment, distance_matrix, matching_breed


# Execute the 'run' function when the script is run directly
if __name__ == "__main__":
    closest_sequence, closest_difference, closest_alignment, distance_matrix, matching_breed = run("../mystery.fa")
    print("Matching Breed:", matching_breed)
    print("Closest DNA Sequence:", closest_sequence)
    print("Closest Alignment:", closest_alignment)
    print("Difference with Closest Sequence:", closest_difference)
    print("Distance Matrix:")
    print(distance_matrix)
