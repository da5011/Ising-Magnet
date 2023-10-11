#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <random>
#include <SDL2/SDL.h>

using namespace std;
typedef vector<vector<int>> lattice; //Defines 2D Vector as a lattice
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<float> dist(0.0, 1.0);

void draw_state(lattice latticeGas, SDL_Renderer* renderer){
    for(int i = 0; i < latticeGas.size(); i++){
        for(int j = 0; j < latticeGas.size(); j++){
            if (latticeGas[i][j] == 0){
                SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
                SDL_RenderDrawPoint(renderer, i, j);
            };
        };
    };
};

void initialize(int lattice_size, lattice &latticeGas){
    latticeGas.resize(lattice_size);
    for (int k=0; k<latticeGas.size(); k++){ 
        latticeGas[k].resize(lattice_size);
    };

    for (int i=0; i<floor(lattice_size*0.5); i++){
        fill(latticeGas[i].begin(), latticeGas[i].end(), 0);
    };

    for (int i=floor(lattice_size*0.5); i<latticeGas.size(); i++){
        fill(latticeGas[i].begin(), latticeGas[i].end(), 1);
    };
};

void valueSwap(lattice &latticeGas, int x1, int y1, int x2, int y2){
    int temp = latticeGas[x1][y1];
    latticeGas[x1][y1] = latticeGas[x2][y2];
    latticeGas[x2][y2] = temp;
}

int main(int argc, char *argv[]){
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    SDL_CreateWindowAndRenderer(100*10, 100*10, 0, &window, &renderer);
    SDL_RenderSetScale(renderer, 10, 10);

    const int lattice_size = 100;
    const int number_of_sweeps = 10000;
    constexpr int swaps_per_sweep = lattice_size*lattice_size;
    lattice latticeGas;

    initialize(lattice_size, latticeGas);
    
    for (int i = 0; i < number_of_sweeps; i++){
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        draw_state(latticeGas, renderer);
        SDL_RenderPresent(renderer);
        SDL_Delay(2);

        for (int j = 0; j < swaps_per_sweep; j++){
            int rand_x = floor(lattice_size*dist(mt));
            int rand_y = floor(lattice_size*dist(mt));
            int neighboor = floor(dist(mt)*4);

            if ((rand_x - 1 > 0) && (neighboor == 0)){
                valueSwap(latticeGas, rand_x, rand_y, rand_x-1, rand_y);
                continue;

            }else if ((rand_y - 1 > 0) && (neighboor == 1)){
                valueSwap(latticeGas, rand_x, rand_y, rand_x, rand_y-1);
                continue;

            }else if ((rand_x + 1 < lattice_size) && (neighboor == 2)){
                valueSwap(latticeGas, rand_x, rand_y, rand_x+1, rand_y);
                continue;

            }else if ((rand_y + 1 < lattice_size) && (neighboor == 3)){
                valueSwap(latticeGas, rand_x, rand_y, rand_x, rand_y+1);
                continue;
            };
        };
    };

    for (int i=0; i<lattice_size; i++){
        for (int j=0; j<lattice_size; j++){
            cout << latticeGas[i][j] << " ";
        }
        cout << "\n";
    }

    return 0;
};
