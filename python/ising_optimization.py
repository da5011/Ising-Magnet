import os
import pandas as pd
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from math import exp
from PIL import Image
import time

######################
#FUNCTION DEFINITIONS#
######################

def startProgram():
    global filePath, simulation_size, temperature, inv_temperature, number_of_sweeps, steps_per_sweep, saveFrames, per_particle
    
    print("**********************")
    print("ISING MAGNET SIMULATOR")
    print("**********************")

    simulation_size = 250 #int(input("Please input the size of the simulation (integer): "))
    temperature = 0.1 #float(input("Please input the temperature of the simulation [0-10] (float): "))
    inv_temperature = 1/temperature
    number_of_sweeps = 1000 #int(input("Please input the number of sweeps for the simulation (integer): "))
    steps_per_sweep = simulation_size*simulation_size
    saveFramesInput = "N" #str(input("Would you like to save simulation snapshots as images? This incurs a performance cost. (Y/N) ")).upper()

    saveFrames = True if saveFramesInput == "Y" else False

    per_particle = 1/steps_per_sweep 

    os.chdir(os.path.dirname(os.path.abspath(__file__))) #Sets CWD to script
    filePath = os.getcwd() + f"\\SimulationData\\{simulation_size}-{temperature}-{number_of_sweeps}"
    if (os.path.exists(filePath + "\\frames") == False):
        os.makedirs(filePath + "\\frames")
    
def initialize(gridSize):
    global energy_array, energyAverage_array, magnitizationSquared_array, \
    energySquared_array, energySquared_average_array, energyAverage_squared_array, \
    heatCapacity_array

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
    
    energy_array[0] = calculateEnergy(isingModel)
    energySquared_array[0] = (energy_array[0])**2
    magnitizationSquared_array[0] = (np.sum(isingModel))**2

    if saveFrames:
        plotModel(isingModel, "Sweep0")
    
    return isingModel

def arrayRunningAverage(array, index):
    return np.average(array[:index+1])

def calculateEnergy(isingModel):
    energy = 0

    for i in range(simulation_size): #Energy calculation w/ PBC
        for j in range(simulation_size):
            currentSpin = isingModel[i,j]
            energy += currentSpin*(isingModel[i,(j+1)%simulation_size] + isingModel[i,(j-1)%simulation_size] + \
                                   isingModel[(i+1)%simulation_size,j] + isingModel[(i-1)%simulation_size,j]) 
    return int(-0.5*energy)

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

def colorSwap(color, isingModel, sweep):
    energySum = 0
    color_dE = 2*isingModel[color]*((np.roll(isingModel,1, axis=0))[color] + (np.roll(isingModel,-1, axis=0))[color] +
                                   (np.roll(isingModel,1, axis=1))[color] + (np.roll(isingModel,-1, axis=1))[color])
    
    color_prob = [True if (exp(-inv_temperature*index) > rand()) else False for index in color_dE]

    isingModel[color] = [-i if j == True else i for i,j in zip(isingModel[color], color_prob)]

    energySum = sum([i if j == True else 0 for i,j in zip(color_dE, color_prob)])

    energy_array[sweep] += energySum

def MonteCarloLoop(number_of_sweeps, isingModel) :
    print("SIMULATION START")

    #Defining grid "colors"
    x,y = np.indices((simulation_size, simulation_size))
    red = np.logical_and(((x+y)%2) == 0, 
                        np.logical_and(x < (simulation_size-1), y < (simulation_size-1)))
    red[simulation_size-1, simulation_size-1] = True

    blue = np.logical_and(((x+y)%2) == 1,
                        np.logical_and(x < (simulation_size-1), y < (simulation_size-1)))

    green = np.logical_and(((x+y)%2) == 0, 
                        np.logical_or(x == (simulation_size-1), y == (simulation_size-1)))
    green[(simulation_size-1),(simulation_size-1)] = False

    yellow = np.logical_and(((x+y)%2) == 1,
                        np.logical_or(x == (simulation_size-1), y == (simulation_size-1)))
    
    for sweep in range(1, number_of_sweeps+1):
        print(f"Sweep{sweep}", end="\r")
        energy_array[sweep] = energy_array[sweep-1]
        colorSwap(red, isingModel, sweep)
        colorSwap(blue, isingModel, sweep)
        colorSwap(green, isingModel, sweep)
        colorSwap(yellow, isingModel, sweep)

def DataAnalysis():
    EnergyAnalysis()
    simulationData = pd.DataFrame({'E':energy_array, '<E>':energyAverage_array, "<E^2>":energySquared_average_array,
                                   "<E>^2":energyAverage_squared_array, "HeatCapacity":heatCapacity_array}) #"<M>^2":magnitizationSquared_array}
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
tic = time.time()
ising = initialize(simulation_size)
MonteCarloLoop(number_of_sweeps, ising)
toc = time.time()
DataAnalysis()


print("Simulation Completed in: " + str(toc-tic) + " seconds.")