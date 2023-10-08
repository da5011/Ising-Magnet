import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from math import floor, exp
from PIL import Image
import os

simulation_size = 128
coupling_constant = 1
temperature = 1
inv_temperature = 1/temperature
number_of_sweeps = 50

saveFrames = True

os.chdir(os.path.dirname(os.path.abspath(__file__))) #Sets CWD to script

filePath = os.getcwd() + "\SimulationData\Temperature=" + str(temperature)
if (os.path.exists(os.getcwd() + "\SimulationData\Temperature=" + str(temperature) + "\\frames") == False):
    os.makedirs(os.getcwd() + "\SimulationData\Temperature=" + str(temperature) + "\\frames")

def plusminusone():
    if rand() > 0.5:
        return 1
    else:
        return -1

def initialize(gridSize):
    print("Initializing Simulation")
    isingModel = np.zeros((gridSize,gridSize))

    with np.nditer(isingModel, op_flags=['readwrite']) as isingIterator:
        for x in isingIterator:
            x[...] = plusminusone()
    
    if saveFrames:
        plotModel(isingModel, "Sweep0")
    return isingModel

def calculateEnergy(isingModel):
    energy = 0
    nrows = len(isingModel)
    ncols = len(isingModel[0])

    for i in range(nrows):
        for j in range(ncols):
            currentSpin = isingModel[i,j]
            energy += currentSpin*(isingModel[i,(j+1)%ncols] + isingModel[i,(j-1)%ncols] + \
                                   isingModel[(i+1)%nrows,j] + isingModel[(i-1)%nrows,j]) #Energy calculation w/ PBC
    return int((-energy*coupling_constant)/2)

def spinFlipEnergyChange(isingModel, x_coord, y_coord):
    nrows = len(isingModel)
    ncols = len(isingModel[0])
    currentSpin = isingModel[x_coord,y_coord]
    deltaE = (2*currentSpin)*(isingModel[x_coord,(y_coord+1)%ncols] + isingModel[x_coord,(y_coord-1)%ncols] + \
                              isingModel[(x_coord+1)%nrows,y_coord] + isingModel[(x_coord-1)%nrows,y_coord]) #Energy calculation w/ PBC
    
    return deltaE

def plotModel(model, fileName):
    plt.imshow(model, cmap= "jet")
    plt.title("2-D Ising Model")
    plt.axis('off')
    plt.savefig(filePath + "\\frames" + fileName, dpi = 800, bbox_inches = "tight")
    plt.close()

def MonteCarloLoop(isingModel, number_of_sweeps):
    steps_per_sweep = len(isingModel)*len(isingModel[0])

    for sweep in range(0, number_of_sweeps):
        print("Computing Sweep #" + str(sweep+1), end="\r")
        for step in range(0, steps_per_sweep):
            rand_x = floor(rand()*len(isingModel))
            rand_y = floor(rand()*len(isingModel[0]))

            if(exp(-inv_temperature*spinFlipEnergyChange(isingModel, rand_x, rand_y)) > rand()):
                isingModel[rand_x, rand_y] = -isingModel[rand_x, rand_y]

        if saveFrames:
            plotModel(isingModel, f"Sweep{sweep+1}")
        print("Simulation Complete")
    return isingModel

my_model = initialize(simulation_size)

#MonteCarloLoop(my_model, number_of_sweeps)