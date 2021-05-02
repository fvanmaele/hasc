#!/usr/bin/env python3
# import matplotlib.pyplot as plt
# import numpy as np
# import csv
import pandas as pd

df = pd.read_csv('transpose_gcc_O3.txt')
fig = df.plot(x='N', y=[1:]).get_figure()

fig.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
fig.savefig("transpose_plot.png", bbox_inches="tight")
