# Import necessary modules from Biopython
from Bio import SeqIO  # Used to parse sequences from files
from Bio.Phylo.TreeConstruction import DistanceCalculator, \
    DistanceTreeConstructor  # Used to calculate distances between sequences and construct trees
from Bio.Phylo import draw  # Used to draw phylogenetic trees
from Bio.SeqRecord import SeqRecord  # Used to represent individual sequences as SeqRecord objects
from Bio.Align import MultipleSeqAlignment  # Used to create multiple sequence alignments
from Bio import pairwise2


def read_dna_files(filepath):
    lines_where_breed_are_present = list()
    dog_breed_with_dna = list()

    # find which lines contain ">gb"
    with open(filepath) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i].strip()
            if line.startswith(">gb"):
                lines_where_breed_are_present.append(i)

    if len(lines_where_breed_are_present) == 1:
        return "".join(lines[1:]).replace("\n", "")
    else:
        # finding out the breed name and the dna sequence
        for i in range(1, len(lines_where_breed_are_present)):
            prev = lines_where_breed_are_present[i - 1]
            curr = lines_where_breed_are_present[i]
            line = lines[prev]
            details = line.split("[")
            for detail in details:
                if detail.startswith("breed"):
                    breed = detail.replace("breed=", "").replace("]", "")
                    dog_breed_with_dna.append([breed, "".join(lines[prev + 1:curr]).replace("\n", "")])

    return dog_breed_with_dna

def compare_alignments(mystery_breed_dna_sequence, known_breed_dna_sequence):
    alignments = pairwise2.align.globalxx(known_breed_dna_sequence, mystery_breed_dna_sequence)
    return alignments

# Define a function to calculate the similarity between a mystery DNA sequence and known dog breeds
def calculate(dog_breed_with_dna, mystery_dna_sequence):
    matching_breed = None
    max_score = 0
    most_compatible_sequence = None
    text = ''
    i = 0
    for breed, dna_sequence in dog_breed_with_dna:
        i += 1
        print(i)
        alignments = compare_alignments(dna_sequence, mystery_dna_sequence)
        for detail in alignments:
            score = detail[2]
            if score > max_score:
                max_score = score
                matching_breed = breed
                most_compatible_sequence = detail
        text += matching_breed + '\n' + str(max_score) + '\n' + str(most_compatible_sequence) + '\n\n'

    with open("results.txt", "w") as f:
        f.write(text)

    print(matching_breed)
    print(max_score)
    print(most_compatible_sequence)
    return matching_breed

# Define a function to read DNA sequences from a file
def read_dna_files_for_tree(filepath):
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


if __name__ == "__main__":
    # read the known dog breeds dna file
    dog_breed_with_dna = read_dna_files("../dog_breeds.fa")

    # read the mystery dog breed dna file
    mystery_dna_sequence = read_dna_files("../mystery.fa")

    # get details for comparison for mystery and known dog dna sequences
    matching_breed = calculate(dog_breed_with_dna, mystery_dna_sequence)

    dog_breeds = read_dna_files_for_tree("../dog_breeds.fa")
    generate_tree(dog_breeds)
