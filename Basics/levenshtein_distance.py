def levenshtein_distance(seq1, seq2):
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    
    matrix = [[0 for _ in range(size_y)] for _ in range(size_x)]
    
    for x in range(size_x):
        matrix[x][0] = x
    for y in range(size_y):
        matrix[0][y] = y
        
    for x in range(1, size_x):
        for y in range(1, size_y):
            if seq1[x-1] == seq2[y-1]:
                matrix[x][y] = matrix[x-1][y-1]
            else:
                matrix[x][y] = min(
                    matrix[x-1][y] + 1,    # Deletion
                    matrix[x][y-1] + 1,    # Insertion
                    matrix[x-1][y-1] + 1   # Substitution
                )
                
    return matrix[size_x - 1][size_y - 1]

if __name__ == "__main__":
    sequence_a = "GATCGGCAT"
    sequence_b = "CAACGGCAG"
    
    distance = levenshtein_distance(sequence_a, sequence_b)
    
    print(f"Sequence 1: {sequence_a}")
    print(f"Sequence 2: {sequence_b}")
    print(f"Levenshtein Distance: {distance}")

    max_len = max(len(sequence_a), len(sequence_b))
    similarity = (1 - distance / max_len) * 100
    print(f"Similarity: {similarity:.2f}%")
