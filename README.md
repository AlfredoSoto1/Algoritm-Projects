# Needleman-Wunsch Algorithm Project

## Overview

This project implements the **Needleman-Wunsch algorithm**, a classic method for global sequence alignment. It is widely used in bioinformatics to align protein or nucleotide sequences, scoring matches, mismatches, and gaps based on a predefined scoring system.

The program reads a CSV file containing pairs of sequences, aligns them using the Needleman-Wunsch algorithm, and outputs the aligned sequences along with their alignment scores.

---

## Features

1. **Global Sequence Alignment**:
   - Computes the optimal alignment of two sequences based on a scoring matrix and gap penalties.
   - Outputs aligned sequences and the alignment score.

2. **Customizable Scoring**:
   - Match score: `+1`.
   - Mismatch penalty: `-1`.
   - Gap penalty: `-2`.

3. **CSV File Integration**:
   - Reads pairs of sequences from a CSV file for alignment.

4. **Alignment Output**:
   - Outputs aligned sequences and their alignment scores to the console.

---

## Files in the Project

1. **`needleman_wunsch.py`**: The primary Python script containing the Needleman-Wunsch implementation.
2. **Input file (CSV)**: A file containing pairs of sequences for alignment (e.g., `sequences.csv`).

---

## Installation and Usage

### Prerequisites

- Python 3.7 or higher.
- The `csv` module (built into Python).
- A terminal or command-line interface.

### Running the Project

1. **Prepare the input file**:
   - Create a CSV file where each row contains a pair of sequences to be aligned.
   - Example format:
     ```
     Sequence1,Sequence2
     AGCT,AGTT
     GAT,CTT
     ```

2. **Run the script**:
   ```bash
   python needleman_wunsch.py <input_file>
