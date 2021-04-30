#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import csv

plt.figure(figsize=(10, 6), dpi=100)
with open('pointer_chasing_results_clang_O1_native.txt') as csvfile:
    reader = csv.reader(csvfile, skipinitialspace=True)
    header = np.array(next(reader))
    max_row_length = len(header) - 1
    
    for index, row in enumerate(reader):
        # Add missing values as NaN
        # https://matplotlib.org/devdocs/gallery/lines_bars_and_markers/masked_demo.html
        row = np.append(row, np.zeros(max_row_length - len(row) + 1) + np.nan)

        # Convert rows to the correct data type
        x = header[1:].astype(float)
        y = row[1:].astype(float)
        plt.loglog(x, y, marker='*', label=row[0])

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig("pointer_chasing.png", bbox_inches="tight")