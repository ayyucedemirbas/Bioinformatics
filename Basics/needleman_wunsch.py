import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-1):
    m, n = len(seq1), len(seq2)
    
    matrix = np.zeros((m + 1, n + 1))
    
    for i in range(m + 1):
        matrix[i][0] = i * gap
    for j in range(n + 1):
        matrix[0][j] = j * gap
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                diagonal = matrix[i-1][j-1] + match
            else:
                diagonal = matrix[i-1][j-1] + mismatch
                
            up = matrix[i-1][j] + gap
            left = matrix[i][j-1] + gap
            
            matrix[i][j] = max(diagonal, up, left)
    
    align1, align2 = "", ""
    i, j = m, n
    
    while i > 0 or j > 0:
        current_score = matrix[i][j]
        
        if i > 0 and j > 0:
            if seq1[i-1] == seq2[j-1]:
                diagonal_score = matrix[i-1][j-1] + match
            else:
                diagonal_score = matrix[i-1][j-1] + mismatch
        else:
            diagonal_score = float('-inf')
            
        up_score = matrix[i-1][j] + gap if i > 0 else float('-inf')
        left_score = matrix[i][j-1] + gap if j > 0 else float('-inf')
        
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
    
    return matrix[m][n], (align1, align2), matrix

def visualize_matrix(matrix, seq1, seq2, title="Needleman-Wunsch Scoring Matrix"):
    fig, ax = plt.subplots(figsize=(max(8, len(seq2)), max(6, len(seq1))))
    
    row_labels = ["-"] + list(seq1)
    col_labels = ["-"] + list(seq2)
    
    sns.heatmap(matrix, annot=True, fmt='.0f', cmap='RdYlBu_r',
                xticklabels=col_labels, yticklabels=row_labels,
                cbar_kws={'label': 'Score'}, ax=ax)
    
    ax.set_title(title, fontsize=14, pad=20)
    ax.set_xlabel('Sequence 2', fontsize=12)
    ax.set_ylabel('Sequence 1', fontsize=12)
    
    plt.tight_layout()
    plt.show()

def print_alignment(align1, align2, score):
    print(f"\nOptimal Alignment (Score: {score}):")
    print("-" * 40)
    print(f"Seq1: {align1}")
    print(f"Seq2: {align2}")
    
    match_line = ""
    for i in range(len(align1)):
        if align1[i] == align2[i] and align1[i] != '-':
            match_line += "|"
        elif align1[i] == '-' or align2[i] == '-':
            match_line += " "
        else:
            match_line += "x"
    print(f"      {match_line}")

if __name__ == "__main__":
    seq1 = "ACGTAG"
    seq2 = "AGTACG"
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Parameters: Match=+2, Mismatch=-1, Gap=-1")

    score, alignment, matrix = needleman_wunsch(seq1, seq2)
    
    print_alignment(alignment[0], alignment[1], score)
    
    visualize_matrix(matrix, seq1, seq2)
  
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
