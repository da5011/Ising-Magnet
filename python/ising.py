import os
import pandas as pd
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from math import floor, exp
from PIL import Image

######################
#FUNCTION DEFINITIONS#
######################

def startProgram():
    global filePath, simulation_size, temperature, inv_temperature, number_of_sweeps, saveFrames
    
    print("**********************")
    print("ISING MAGNET SIMULATOR")
    print("**********************")

    simulation_size = int(input("Please input the size of the simulation (integer): "))
    temperature = float(input("Please input the temperature of the simulation [0-10] (float): "))
    inv_temperature = 1/temperature
    number_of_sweeps = int(input("Please input the number of sweeps for the simulation (integer): "))
    saveFramesInput = str(input("Would you like to save simulation snapshots as images? This incurs a performance cost. (Y/N) ")).upper()

    saveFrames = True if saveFramesInput == "Y" else False

    os.chdir(os.path.dirname(os.path.abspath(__file__))) #Sets CWD to script
    filePath = os.getcwd() + f"\SimulationData\{simulation_size}-{temperature}-{number_of_sweeps}"
    if (os.path.exists(filePath + "\\frames") == False):
        os.makedirs(filePath + "\\frames")

def plusminusone():
    if rand() > 0.5:
        return 1
    else:
        return -1

def initialize(gridSize):
    global EnergyPerSweep

    print("Initializing Simulation")

    isingModel = np.zeros((gridSize,gridSize))
    with np.nditer(isingModel, op_flags=['readwrite']) as isingIterator:
        for x in isingIterator:
            x[...] = plusminusone()

    EnergyPerSweep = np.zeros(number_of_sweeps+1)
    EnergyPerSweep[0] = calculateEnergy(isingModel)

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
    return int(-energy/2)

def spinFlipEnergyChange(isingModel, x_coord, y_coord):
    nrows = len(isingModel)
    ncols = len(isingModel[0])
    currentSpin = isingModel[x_coord,y_coord]
    deltaE = (2*currentSpin)*(isingModel[x_coord,(y_coord+1)%ncols] + isingModel[x_coord,(y_coord-1)%ncols] + \
                              isingModel[(x_coord+1)%nrows,y_coord] + isingModel[(x_coord-1)%nrows,y_coord]) #Energy calculation w/ PBC
    
    return deltaE

def plotModel(model, fileName):
    plt.imshow(model, cmap= "jet")
    plt.axis('off')
    plt.savefig(filePath + "\\frames\\" + fileName + ".png", dpi = 800, bbox_inches = "tight")
    plt.close()

def MonteCarloLoop(isingModel, number_of_sweeps):
    print("SIMULATION START")
    steps_per_sweep = len(isingModel)*len(isingModel[0])

    for sweep in range(0, number_of_sweeps):
        print("Computing Sweep #" + str(sweep+1), end="\r")
        EnergyPerSweep[sweep+1] += EnergyPerSweep[sweep]

        for step in range(0, steps_per_sweep):
            rand_x = floor(rand()*len(isingModel))
            rand_y = floor(rand()*len(isingModel[0]))
            spinFlipCandidate = spinFlipEnergyChange(isingModel, rand_x, rand_y)

            if(exp(-inv_temperature*spinFlipCandidate) > rand()):
                isingModel[rand_x, rand_y] = -isingModel[rand_x, rand_y]
                EnergyPerSweep[sweep+1] += spinFlipCandidate

        if saveFrames:
            plotModel(isingModel, f"Sweep{sweep+1}")
    print("Simulation Complete")
    return isingModel

def DataAnalysis(energyArray):
    simulationData = pd.DataFrame(energyArray, columns=["Energy"])
    simulationData.to_excel(filePath + "\SimulationResults.xlsx")

    if saveFrames:
        print("Generating gif...")
        frames = []
        for i in range(0, number_of_sweeps+1):
            frames.append(Image.open(filePath + "\\frames\\" + f"Sweep{i}.png"))
        animation = frames[0]
        animation.save(filePath + "\simulation.gif", format= "gif", append_images= frames, save_all= True, optimize= True, loop= 0)

    print(f"Analysis Complete. Data has been saved to {filePath}")

######################
####PROGRAM START#####
######################

startProgram()
my_model = initialize(simulation_size)
MonteCarloLoop(my_model, number_of_sweeps)
DataAnalysis(EnergyPerSweep)