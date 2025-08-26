import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-1):
    m, n = len(seq1), len(seq2)
    
    matrix = np.zeros((m + 1, n + 1))
    
    max_score = 0
    max_pos = (0, 0)
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                diagonal = matrix[i-1][j-1] + match
            else:
                diagonal = matrix[i-1][j-1] + mismatch
                
            up = matrix[i-1][j] + gap
            left = matrix[i][j-1] + gap
            
            matrix[i][j] = max(0, diagonal, up, left)
            
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)
    
    align1, align2 = "", ""
    i, j = max_pos
    
    while i > 0 and j > 0 and matrix[i][j] > 0:
        current_score = matrix[i][j]
        
        if seq1[i-1] == seq2[j-1]:
            diagonal_score = matrix[i-1][j-1] + match
        else:
            diagonal_score = matrix[i-1][j-1] + mismatch
            
        up_score = matrix[i-1][j] + gap if i > 0 else 0
        left_score = matrix[i][j-1] + gap if j > 0 else 0
        
        if current_score == diagonal_score:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif current_score == up_score:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
    
    start_pos = (i, j)
    return max_score, (align1, align2), matrix, start_pos

def visualize_matrix(matrix, seq1, seq2, max_pos, title="Smith-Waterman Scoring Matrix"):
    fig, ax = plt.subplots(figsize=(max(8, len(seq2)), max(6, len(seq1))))
    
    row_labels = ["-"] + list(seq1)
    col_labels = ["-"] + list(seq2)

    sns.heatmap(matrix, annot=True, fmt='.0f', cmap='RdYlBu_r',
                xticklabels=col_labels, yticklabels=row_labels,
                cbar_kws={'label': 'Score'}, ax=ax)

    ax.add_patch(plt.Rectangle((max_pos[1]-0.5, max_pos[0]-0.5), 1, 1, 
                              fill=False, edgecolor='red', lw=3))
    
    ax.set_title(title, fontsize=14, pad=20)
    ax.set_xlabel('Sequence 2', fontsize=12)
    ax.set_ylabel('Sequence 1', fontsize=12)
    
    plt.tight_layout()
    plt.show()

def print_alignment(align1, align2, score, seq1, seq2, start_pos):
    print(f"\nOptimal Local Alignment (Score: {score}):")
    print("-" * 50)
    
    end_i, end_j = start_pos[0] + len(align1.replace('-', '')), start_pos[1] + len(align2.replace('-', ''))
    
    print(f"Seq1 [{start_pos[0]}:{end_i}]: {align1}")
    print(f"Seq2 [{start_pos[1]}:{end_j}]: {align2}")
    
    match_line = ""
    for i in range(len(align1)):
        if align1[i] == align2[i] and align1[i] != '-':
            match_line += "|"
        elif align1[i] == '-' or align2[i] == '-':
            match_line += " "
        else:
            match_line += "x"
    print(f"{'':17} {match_line}")

def find_all_local_alignments(matrix, seq1, seq2, threshold=0.8):
    max_score = np.max(matrix)
    threshold_score = max_score * threshold
    
    alignments = []
    positions = np.where(matrix >= threshold_score)
    
    for i, j in zip(positions[0], positions[1]):
        if matrix[i][j] >= threshold_score:
            alignments.append((matrix[i][j], (i, j)))
    
    alignments.sort(reverse=True, key=lambda x: x[0])
    return alignments[:5]

if __name__ == "__main__":
    seq1 = "ACGTACGTACGT"
    seq2 = "TGCACGTACGAA"
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Parameters: Match=+2, Mismatch=-1, Gap=-1")
    
    score, alignment, matrix, start_pos = smith_waterman(seq1, seq2)  
    print_alignment(alignment[0], alignment[1], score, seq1, seq2, start_pos)
    
    visualize_matrix(matrix, seq1, seq2, (score, matrix).index(np.max(matrix)) 
                    if hasattr(matrix, 'index') else np.unravel_index(np.argmax(matrix), matrix.shape))
    
    top_alignments = find_all_local_alignments(matrix, seq1, seq2)
    print(f"\nTop local alignments:")
    for i, (score, pos) in enumerate(top_alignments[:3], 1):
        print(f"{i}. Score: {score:.0f} at position {pos}")
    
    print(f"\nScoring Matrix:")
    print("    ", end="")
    for char in "-" + seq2:
        print(f"{char:4}", end="")
    print()
    
    for i, row in enumerate(matrix):
        if i == 0:
            print(f"-   ", end="")
        else:
            print(f"{seq1[i-1]}   ", end="")
        for val in row:
            print(f"{int(val):4}", end="")
        print()
