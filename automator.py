import os
import numpy as np

START = 0.0001
END = 0.01

## WORKS ON LINUX ONLY
def execute_simulation(oil_thickeness):
    print(f'Running simulation for oil thickness = {oil_thickeness:.4f}')
    os.system(f'g++ -O2 --std=c++17 -o output ray_tracer_final.cpp && ./output {oil_thickeness:.4f} data/data_{oil_thickeness:.4f}.txt')

def run_automator():
    for oil_thickeness in np.linspace(START, END, num=30, endpoint=True):
        execute_simulation(oil_thickeness)


if __name__ == '__main__':
    run_automator()