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
    global filePath, simulation_size, temperature, inv_temperature, number_of_sweeps, steps_per_sweep, saveFrames, per_particle
    
    print("**********************")
    print("ISING MAGNET SIMULATOR")
    print("**********************")

    simulation_size = int(input("Please input the size of the simulation (integer): "))
    temperature = float(input("Please input the temperature of the simulation [0-10] (float): "))
    inv_temperature = 1/temperature
    number_of_sweeps = int(input("Please input the number of sweeps for the simulation (integer): "))
    steps_per_sweep = simulation_size*simulation_size
    saveFramesInput = str(input("Would you like to save simulation snapshots as images? This incurs a performance cost. (Y/N) ")).upper()

    saveFrames = True if saveFramesInput == "Y" else False

    per_particle = 1/steps_per_sweep 

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

def calculateEnergy(isingModel): #heat capacity is the variance in the energy w T
    energy = 0
    nrows = len(isingModel)
    ncols = len(isingModel[0])

    for i in range(nrows):
        for j in range(ncols):
            currentSpin = isingModel[i,j]
            energy += currentSpin*(isingModel[i,(j+1)%ncols] + isingModel[i,(j-1)%ncols] + \
                                   isingModel[(i+1)%nrows,j] + isingModel[(i-1)%nrows,j]) #Energy calculation w/ PBC
    return int(-0.5*energy)

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

def compute_AVG_E(energyArray):
    AVG_e = np.zeros(len(energyArray))
    for i in range(0, len(energyArray)):
        AVG_e[i] = np.mean(energyArray[:i+1])
    return AVG_e

def compute_E_minus_AVG_E_SQUARED(energyArray, avgEnergyArray):
    E_minus_AVG_E_SQUARED_Array = np.zeros(len(energyArray))
    for i in range(0, len(energyArray)):
        E_minus_AVG_E_SQUARED_Array[i] = (energyArray[i]-avgEnergyArray[i])**2
    return E_minus_AVG_E_SQUARED_Array

def compute_heatCapacity(bigAssArray):
    return (inv_temperature**2)*compute_AVG_E(bigAssArray)*per_particle

# def compute_E2_AVG(energyArray):
#     E2_AVG_Array = np.zeros(len(energyArray))
#     E2_Array = np.zeros(len(energyArray)) test

#     for i in range(0, len(E2_AVG_Array)):   
#         E2_Array[i] = energyArray[i]*energyArray[i]

#     E2_AVG_Array = compute_AVG_E(E2_Array)
#     return E2_AVG_Array

# def compute_AVG_E2(averageEnergyArray):
#     AVG_E2_Array = np.zeros(len(averageEnergyArray))

#     for i in range(0, len(averageEnergyArray)): 
#         AVG_E2_Array[i] = averageEnergyArray[i]*averageEnergyArray[i]
#     return AVG_E2_Array


def MonteCarloLoop(isingModel, number_of_sweeps):
    print("SIMULATION START")

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

def DataAnalysis():
    averageEnergy = compute_AVG_E(EnergyPerSweep)
    energyDifference = compute_E_minus_AVG_E_SQUARED(EnergyPerSweep, averageEnergy)
    heatCapacity = compute_heatCapacity(energyDifference)
    simulationData = pd.DataFrame({"Energy":EnergyPerSweep, "<E>":averageEnergy,
                                   "[E-<E>]^2":energyDifference, "Heat Capacity":heatCapacity})
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
DataAnalysis()