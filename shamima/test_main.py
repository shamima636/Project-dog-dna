from .main import *

dog_breed_with_dna = read_dna_files("../dog_breeds.fa")
mystery_dna_sequence = read_dna_files("../mystery_test.fa")
matching_breed = calculate(dog_breed_with_dna, mystery_dna_sequence)
if matching_breed == "Portuguese Warren dog, small size, smooth hair":
    print("Testing successful")
else:
    print("Testing failed.")