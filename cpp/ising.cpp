#include <iostream>
#include <math.h>
#include <random>
#include <SDL2/SDL.h>

using namespace std;
typedef vector<vector<int>> lattice; //Defines 2D Vector as a lattice
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<float> dist(0.0, 1.0);

int plusminusone();
int superMod(int x, int y);
void initializeSimulation(int latticeSize, lattice &latticeModel);
void drawFrame(lattice latticeModel, SDL_Renderer* renderer);
int calculateEnergy(lattice latticeModel);
int spinFlipEnergyChange(lattice latticeModel, int x_coord, int y_coord);

int main(int argc, char *argv[]){
    const int latticeSize = 100;
    const int number_of_sweeps = 1000;
    constexpr int steps_per_sweep = latticeSize*latticeSize;
    float temperature = 2.3;
    float inv_temperature = 1/temperature;
    lattice isingMagnet;
    
    const int windowMultiplier = 10;
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    SDL_CreateWindowAndRenderer(latticeSize*windowMultiplier, latticeSize*windowMultiplier, 0, &window, &renderer);
    SDL_RenderSetScale(renderer, windowMultiplier, windowMultiplier);
    SDL_SetWindowTitle(window, "Ising Magnet Monte Carlo");

    initializeSimulation(latticeSize, isingMagnet);

    for (int n=0; n < number_of_sweeps; n++){
        drawFrame(isingMagnet, renderer);
        for (int m=0; m < steps_per_sweep; m++){
            int rand_x = floor(latticeSize*dist(rd));
            int rand_y = floor(latticeSize*dist(rd));

            int spinFlipCandidate = spinFlipEnergyChange(isingMagnet, rand_x, rand_y);

            if(exp(-inv_temperature*spinFlipCandidate) > dist(rd)){
                isingMagnet[rand_x][rand_y] = -isingMagnet[rand_x][rand_y];
            }
        }
    }

    for (int i=0; i<latticeSize; i++){
        for (int j=0; j<latticeSize; j++){
            cout << isingMagnet[i][j] << " ";
        }
        cout << "\n";
    }
    return 0;
};

void initializeSimulation(int latticeSize, lattice &latticeModel){
    latticeModel.resize(latticeSize);
    for (int i = 0; i < latticeSize; i++){
        latticeModel[i].resize(latticeSize);
    };

    for (int i = 0; i < latticeSize; i++){
        for(int j = 0; j < latticeSize; j++){
            latticeModel[i][j] = plusminusone();
        };
    };
    int intial_energy = calculateEnergy(latticeModel);
}

void drawFrame(lattice latticeModel, SDL_Renderer* renderer){
    SDL_SetRenderDrawColor(renderer, 10, 45, 70, 255);
    SDL_RenderClear(renderer);

    for(int i = 0; i < latticeModel.size(); i++){
        for(int j = 0; j < latticeModel.size(); j++){
            if (latticeModel[i][j] == -1){
                SDL_SetRenderDrawColor(renderer, 151, 209, 187, 255);
                SDL_RenderDrawPoint(renderer, i, j);
            };
        };
    };

    SDL_RenderPresent(renderer);
}

int calculateEnergy(lattice latticeModel){
    int size = latticeModel.size();
    int energy = 0;
    int sum = 0;
    for (int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            energy += latticeModel[i][j]*(latticeModel[i][superMod(j+1,size)] + latticeModel[i][superMod(j-1,size)] + 
                                          latticeModel[superMod(i+1,size)][j] + latticeModel[superMod(i-1,size)%size][j]);
        };
    };
    return (-energy/2);
}

int spinFlipEnergyChange(lattice latticeModel, int x_coord, int y_coord){
    int size = latticeModel.size();
    return 2*latticeModel[x_coord][y_coord]*
    (latticeModel[x_coord][superMod(y_coord+1,size)] + latticeModel[x_coord][superMod(y_coord-1,size)] + 
     latticeModel[superMod(x_coord+1,size)][y_coord] + latticeModel[superMod(x_coord-1,size)%size][y_coord]);
}

int plusminusone(){
    return ((dist(rd)) < 0.5) ? -1 : 1;
}

int superMod(int x, int y){
    return ((x %= y) < 0) ? x+y : x;
}