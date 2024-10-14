import sys
import csv

MATCH = 1
MISSMATCH = -1

def scoring_matrix(seq_A, seq_B) -> int:
    return MATCH if seq_A == seq_B else MISSMATCH

def fill_needleman_matrix(gap_penalty, m, n, seq_A, seq_B):
    # Initialize matrix all zeros
    F = []
    for i in range(m + 1):
        F.append([0] * (n + 1))

    # Initialize the edges of the matrix
    for i in range(m + 1):
        F[i][0] = i * gap_penalty
    for j in range(n + 1):
        F[0][j] = j * gap_penalty

    # Fill the scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Find value for scoring matrix
            score = F[i-1][j-1] + scoring_matrix(seq_A[i - 1], seq_B[j - 1])
            
            # Find the delete and instert with the gap penalty
            delete = F[i-1][j] + gap_penalty
            insert = F[i][j-1] + gap_penalty

            # Find the max between the three
            F[i][j] = max(score, delete, insert)

    return F # Return the filled out matrix

def backtrack(m, n, gap_penalty, F, seq_A, seq_B):
    # Traceback to get the alignment
    i = m
    j = n
    aligned_A = ""
    aligned_B = ""

    while i > 0 and j > 0:
        # Backtracking directions
        curr_score = F[i    ][j    ]
        diag_score = F[i - 1][j - 1]
        up_score   = F[i - 1][j    ]
        left_score = F[i    ][j - 1]

        # If the two letters match, mark the score correspondingly
        score = scoring_matrix(seq_A[i-1], seq_B[j-1])

        if curr_score == left_score + gap_penalty:
            # Go left
            aligned_A = "-" + aligned_A
            aligned_B = seq_B[j - 1] + aligned_B
            j -= 1
        elif curr_score == up_score + gap_penalty:
            # Go up
            aligned_A = seq_A[i - 1] + aligned_A
            aligned_B = "-" + aligned_B
            i -= 1
        elif curr_score == diag_score + score:
            # Go diagonaly
            aligned_A = seq_A[i - 1] + aligned_A
            aligned_B = seq_B[j - 1] + aligned_B
            i -= 1
            j -= 1

    # Handle remaining gaps
    while i > 0:
        aligned_A = seq_A[i-1] + aligned_A
        aligned_B = "-" + aligned_B
        i -= 1
    while j > 0:
        aligned_A = "-" + aligned_A
        aligned_B = seq_B[j-1] + aligned_B
        j -= 1

    alignment_score = F[m][n]
    return aligned_A, aligned_B, alignment_score

def needleman_wunsch(seq_A, seq_B):
    # Scoring Matrix
    gap_penalty = -2

    # Dimensions of the matrix
    m = len(seq_A)
    n = len(seq_B)
    
    F = fill_needleman_matrix(gap_penalty, m, n, seq_A, seq_B)

    # for row in F:
    #     line = ""
    #     for value in row:
    #         line += f"{value:4}"
    #     print(line)

    return backtrack(m, n, gap_penalty, F, seq_A, seq_B)

if len(sys.argv) > 1:
    input_file = sys.argv[1]

    with open(input_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader, None) # Skip header 
        for row in reader:
            if len(row) != 2:
                print("Invalid input format.")
                continue
            
            seq_A, seq_B = row
            alignedA, alignedB, score = needleman_wunsch(seq_A, seq_B)

            # Output result for each pair
            print(f"{alignedA} {alignedB} {score}")
