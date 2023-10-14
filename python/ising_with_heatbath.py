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
    filePath = os.getcwd() + f"\\SimulationData\\{simulation_size}-{temperature}-{number_of_sweeps}"
    if (os.path.exists(filePath + "\\frames") == False):
        os.makedirs(filePath + "\\frames")
    
def initialize(gridSize):
    global energy_array, energyAverage_array, magnitizationSquared_array, \
    energySquared_array, energySquared_average_array, energyAverage_squared_array, \
    heatCapacity_array, isingModel

    print("Initializing Simulation")

    isingModel = np.zeros((gridSize,gridSize))
    with np.nditer(isingModel, op_flags=['readwrite']) as isingIterator:
        for x in isingIterator:
            x[...] = 1 if (rand() > 0.5) else -1

    if True:
        energy_array = np.zeros(number_of_sweeps+1)
        energySquared_array = np.zeros(number_of_sweeps+1)
        magnitizationSquared_array = np.zeros(number_of_sweeps+1)
        energySquared_average_array = np.zeros(number_of_sweeps+1)
        energyAverage_squared_array = np.zeros(number_of_sweeps+1)
        heatCapacity_array = np.zeros(number_of_sweeps+1)
        energyAverage_array = np.zeros(number_of_sweeps+1)
    
    
    energy_array[0] = calculateEnergy()
    energySquared_array[0] = (energy_array[0])**2
    magnitizationSquared_array[0] = (np.sum(isingModel))**2

    if saveFrames:
        plotModel(isingModel, "Sweep0")

def arrayRunningAverage(array, index):
    return np.average(array[:index+1])

def calculateEnergy(): #heat capacity is the variance in the energy w T
    energy = 0
    nrows = len(isingModel)
    ncols = len(isingModel[0])

    for i in range(nrows):
        for j in range(ncols):
            currentSpin = isingModel[i,j]
            energy += currentSpin*(isingModel[i,(j+1)%ncols] + isingModel[i,(j-1)%ncols] + \
                                   isingModel[(i+1)%nrows,j] + isingModel[(i-1)%nrows,j]) #Energy calculation w/ PBC
    return int(-0.5*energy)

def spinFlipEnergyChange(x_coord, y_coord):

    currentSpin = isingModel[x_coord,y_coord]
    deltaE = (2*currentSpin)*(isingModel[x_coord,(y_coord+1)%ncols] + isingModel[x_coord,(y_coord-1)%ncols] + \
                              isingModel[(x_coord+1)%nrows,y_coord] + isingModel[(x_coord-1)%nrows,y_coord]) #Energy calculation w/ PBC
    
    return deltaE

def plotModel(model, fileName):
    plt.imshow(model, cmap= "jet")
    plt.axis('off')
    plt.savefig(filePath + "\\frames\\" + fileName + ".png", dpi = 800, bbox_inches = "tight")
    plt.close()

def EnergyAnalysis():
    for i in range(0,number_of_sweeps+1):
        energySquared_array[i] = (energy_array[i])**2
        energyAverage_array[i] = arrayRunningAverage(energy_array, i)
        energySquared_average_array[i] = arrayRunningAverage(energySquared_array, i)
        energyAverage_squared_array[i] = (energyAverage_array[i])**2
        heatCapacity_array[i] = (inv_temperature**2)*(per_particle)*(energySquared_average_array[i]-energyAverage_squared_array[i])

def MonteCarloLoop(number_of_sweeps):
    print("SIMULATION START")
    magnetization = np.sum(isingModel)

    x,y = np.indices((simulation_size,simulation_size))
    red = np.logical_and((x+y) % 2==0   , np.logical_and(x<simulation_size-1, y<simulation_size-1))
    red[simulation_size-1,simulation_size-1] = True

    blue = np.logical_and((x+y) % 2==1  , np.logical_and(x<simulation_size-1, y<simulation_size-1))

    green = np.logical_and((x+y) % 2==0 , np.logical_or(x==simulation_size-1, y==simulation_size-1))
    green[simulation_size-1,simulation_size-1] = False

    yellow = np.logical_and((x+y) % 2==1, np.logical_or(x==simulation_size-1, y==simulation_size-1))

    dE_RED = np.array(len(red))
    dE_RED[i] = (-2*red[i])
    for sweep in range(0, number_of_sweeps):
        print("Computing Sweep #" + str(sweep+1), end="\r")
        energy_array[sweep+1] += energy_array[sweep]

        if saveFrames:
            plotModel(isingModel, f"Sweep{sweep+1}")

    print("\nSimulation Complete")
    return isingModel

def DataAnalysis():
    EnergyAnalysis()
    simulationData = pd.DataFrame({'E':energy_array, '<E>':energyAverage_array, "<E^2>":energySquared_average_array,
                                   "<E>^2":energyAverage_squared_array, "HeatCapacity":heatCapacity_array, "<M>^2":magnitizationSquared_array})
    simulationData.to_excel(filePath + f"\\{simulation_size}-{temperature}-{number_of_sweeps}.xlsx")

    if saveFrames:
        print("Generating gif...")
        frames = []
        for i in range(0, number_of_sweeps+1):
            frames.append(Image.open(filePath + "\\frames\\" + f"Sweep{i}.png"))
        animation = frames[0]
        animation.save(filePath + "\\simulation.gif", format= "gif", append_images= frames, save_all= True, optimize= True, loop= 0)

    print(f"Analysis Complete. Data has been saved to {filePath}")

######################
####PROGRAM START#####
######################

startProgram()
initialize(simulation_size)
MonteCarloLoop(number_of_sweeps)
DataAnalysis()