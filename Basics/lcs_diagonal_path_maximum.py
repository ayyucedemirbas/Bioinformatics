def lcs_max_diagonal_path(seq1, seq2):
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    diag = [[0] * (n + 1) for _ in range(m + 1)]
  
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
                diag[i][j] = diag[i-1][j-1] + 1
            else:
                if dp[i-1][j] > dp[i][j-1]:
                    dp[i][j] = dp[i-1][j]
                    diag[i][j] = diag[i-1][j]
                elif dp[i-1][j] < dp[i][j-1]:
                    dp[i][j] = dp[i][j-1]
                    diag[i][j] = diag[i][j-1]
                else:
                    dp[i][j] = dp[i-1][j]
                    if diag[i-1][j] >= diag[i][j-1]:
                        diag[i][j] = diag[i-1][j]
                    else:
                        diag[i][j] = diag[i][j-1]
    lcs = []
    path = []
    i, j = m, n
    
    while i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1]:
            lcs.append(seq1[i-1])
            path.append(('diagonal', i-1, j-1))
            i -= 1
            j -= 1
        else:
            if dp[i-1][j] > dp[i][j-1]:
                path.append(('up', i-1, j))
                i -= 1
            elif dp[i-1][j] < dp[i][j-1]:
                path.append(('left', i, j-1))
                j -= 1
            else:
                if diag[i-1][j] >= diag[i][j-1]:
                    path.append(('up', i-1, j))
                    i -= 1
                else:
                    path.append(('left', i, j-1))
                    j -= 1
  
    while i > 0:
        path.append(('up', i-1, j))
        i -= 1
    while j > 0:
        path.append(('left', i, j-1))
        j -= 1
    
    lcs.reverse()
    path.reverse()
    
    diagonal_count = sum(1 for move in path if move[0] == 'diagonal')
    
    return ''.join(lcs), diagonal_count, path

def visualize_path(seq1, seq2, path):
    m, n = len(seq1), len(seq2)
    grid = [['.' for _ in range(n + 1)] for _ in range(m + 1)]
    
    for move_type, i, j in path:
        if move_type == 'diagonal':
            grid[i][j] = 'D'
        elif move_type == 'up':
            grid[i][j] = '↑'
        elif move_type == 'left':
            grid[i][j] = '←'
    
    print("    ", end="")
    for j, char in enumerate(" " + seq2):
        print(f"{char:2}", end="")
    print()
    
    for i in range(m + 1):
        if i == 0:
            print("  ", end="")
        else:
            print(f"{seq1[i-1]} ", end="")
        
        for j in range(n + 1):
            print(f"{grid[i][j]:2}", end="")
        print()


def main():
    examples = [
        ("ABCDGH", "AEDFHR"),
        ("AGGTAB", "GXTXAYB"),
        ("ATCGATCG", "ATGCATGC"),
    ]
    
    for seq1, seq2 in examples:
        print(f"Sequences: '{seq1}' and '{seq2}'")
        lcs, diag_count, path = lcs_max_diagonal_path(seq1, seq2)
        
        print(f"LCS: '{lcs}' (length: {len(lcs)})")
        print(f"Diagonal edges: {diag_count}")
        print(f"Total moves: {len(path)}")
        
        diagonal_moves = [p for p in path if p[0] == 'diagonal']
        if diagonal_moves:
            print("Diagonal matches:", end=" ")
            for _, i, j in diagonal_moves:
                print(f"({seq1[i]},{i},{j})", end=" ")
            print()
        
        print("Path visualization:")
        visualize_path(seq1, seq2, path)
        print("-" * 40)


if __name__ == "__main__":
    main()
