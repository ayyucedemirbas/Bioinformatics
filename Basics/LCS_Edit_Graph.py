import matplotlib.pyplot as plt
import numpy as np

s1 = "GATTACA"
s2 = "GCATGCU"
m, n = len(s1), len(s2)

# Compute LCS DP table
dp = np.zeros((m+1, n+1), dtype=int)
for i in range(1, m+1):
    for j in range(1, n+1):
        if s1[i-1] == s2[j-1]:
            dp[i, j] = dp[i-1, j-1] + 1
        else:
            dp[i, j] = max(dp[i-1, j], dp[i, j-1])

path_edges = []
i, j = m, n
while i > 0 or j > 0:
    if i > 0 and j > 0 and s1[i-1] == s2[j-1] and dp[i, j] == dp[i-1, j-1] + 1:
        path_edges.append(((i, j), (i-1, j-1)))
        i -= 1
        j -= 1
    elif i > 0 and dp[i, j] == dp[i-1, j]:
        path_edges.append(((i, j), (i-1, j)))
        i -= 1
    else:
        path_edges.append(((i, j), (i, j-1)))
        j -= 1

path_edges = path_edges

fig, ax = plt.subplots(figsize=(10, 8))

for i in range(m+1):
    for j in range(n+1):
        ax.plot(j, m-i, marker='o', markersize=3)
for i in range(m+1):
    for j in range(n+1):
        if j < n:
            ax.arrow(j, m-i, 0.8, 0, head_width=0.12, head_length=0.2, length_includes_head=True, linewidth=0.5)
        if i < m:
            ax.arrow(j, m-i, 0, -0.8, head_width=0.12, head_length=0.2, length_includes_head=True, linewidth=0.5)

for i in range(1, m+1):
    for j in range(1, n+1):
        if s1[i-1] == s2[j-1]:
            ax.arrow(j-1, m-(i-1), 0.8, -0.8, head_width=0.12, head_length=0.2, length_includes_head=True, linewidth=0.6, linestyle='--')

path_edges_forward = path_edges[::-1]
for (i1,j1),(i2,j2) in path_edges_forward:
    x1, y1 = j1, m - i1
    x2, y2 = j2, m - i2
    ax.plot([x1, x2], [y1, y2], linewidth=3)

coords = [(m - e[0][0], e[0][1]) for e in path_edges_forward] + [(m - 0, 0)]

nodes_xy = [(e[1][1], m - e[1][0]) for e in path_edges_forward]

nodes_xy = [(0, m)] + nodes_xy + [(n, 0)]
xs, ys = zip(*nodes_xy)
ax.plot(xs, ys, marker='s', linestyle='None', markersize=6)

for j, ch in enumerate(s2, start=1):
    ax.text(j, m+0.5, ch, ha='center', va='bottom', fontsize=12)
ax.text(0, m+0.5, 'j', ha='center', va='bottom', fontsize=11)
for i, ch in enumerate(s1, start=1):
    ax.text(-0.5, m-(i-0.5), ch, ha='right', va='center', fontsize=12)
ax.text(-0.5, m+0.1, 'i', ha='right', va='center', fontsize=11)


lcs_len = dp[m, n]
i, j = m, n
lcs_chars = []
while i>0 and j>0:
    if s1[i-1] == s2[j-1]:
        lcs_chars.append(s1[i-1])
        i -= 1; j -= 1
    elif dp[i-1, j] >= dp[i, j-1]:
        i -= 1
    else:
        j -= 1
lcs = ''.join(reversed(lcs_chars))

ax.set_title(f"DNA Sequence Alignment as a Path in the Edit Graph\ns1={s1}, s2={s2}, LCS length={lcs_len}, LCS='{lcs}'")
ax.set_xlim(-1, n+1)
ax.set_ylim(-1, m+1)
ax.set_xticks(range(0, n+1))
ax.set_yticks(range(0, m+1))
ax.set_xlabel("Columns j (sequence s2)")
ax.set_ylabel("Rows i (sequence s1)")
ax.set_aspect('equal')
ax.grid(False)
plt.tight_layout()

out_path = "dna_edit_graph_lcs.png"
plt.savefig(out_path, dpi=200, bbox_inches='tight')
