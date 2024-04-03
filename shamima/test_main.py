from main import *

closest_sequence, closest_difference, closest_alignment, distance_matrix, matching_breed = run("../mystery_test.fa")
if matching_breed == "Portuguese Warren dog, small size, smooth hair":
    print("Testing successful")
else:
    print("Testing failed.")
