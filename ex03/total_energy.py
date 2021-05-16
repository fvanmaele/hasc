#!/usr/bin/env pvpython
import pandas as pd
import numpy as np
from paraview.simple import *

def potential_energy(mass, x0, x1, x2, G=1.0):
    x0_i, x0_j = np.meshgrid(x0, x0)
    x1_i, x1_j = np.meshgrid(x1, x1)
    x2_i, x2_j = np.meshgrid(x2, x2)
    m_i, m_j = np.meshgrid(mass, mass)

    dist = np.sqrt((x0_j - x0_i)**2 
                 + (x1_j - x1_i)**2 
                 + (x2_j - x2_i)**2)
    mask = dist > 0 # dist == 0 iff i == j
    return -1/2 * G * np.sum(m_j[mask] * m_i[mask] / dist[mask])


def kinetic_energy(mass, v0, v1, v2):
    dist_sq = v0**2  + v1**2 + v2**2
    return 1/2 * np.sum(mass * dist_sq)


def main():
    time_steps = range(0, 10) # test_data_points_%d.csv
    G = 1.0 # value for gamma
    E_pot_a = np.zeros(len(time_steps))
    E_kin_a = np.zeros(len(time_steps))
    basename = 'test_threaded' # XXX: take basename from command line

    for t in time_steps:
        basevtk = LegacyVTKReader(FileNames=['vtk/{}_{:06d}.vtk'.format(basename, t)])
        SaveData('csv/{}_{:06d}.csv'.format(basename, t), proxy=basevtk)

        df = pd.read_csv('csv/{}_{:06d}.csv'.format(basename, t))
        mass = np.array(df['mass'])

        x0 = np.array(df['Points:0'])
        x1 = np.array(df['Points:1'])
        x2 = np.array(df['Points:2'])
        #x = np.vstack((x0, x1, x2)).transpose()

        v0 = np.array(df['velocity:0'])
        v1 = np.array(df['velocity:1'])
        v2 = np.array(df['velocity:2'])
        #v = np.vstack((v0, v1, v2)).transpose()

        E_pot_a[t] = potential_energy(mass, x0, x1, x2, G)
        E_kin_a[t] = kinetic_energy(mass, v0, v1, v2)

    print(E_pot_a + E_kin_a)


if __name__ == "__main__":
    main()

