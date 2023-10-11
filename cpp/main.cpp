#include <iostream>
#include <math.h>
#include <random>

using namespace std;
typedef vector<vector<int>> lattice; //Defines 2D Vector as a lattice
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<float> dist(0.0, 1.0);

int plusminusone();
int superMod(int x, int y);
void initializeSimulation(int latticeSize, lattice &latticeModel);
void initializeRender();
int calculateEnergy(lattice latticeModel);
int spinFlipEnergyChange(lattice latticeModel, int x_coord, int y_coord);

int main(){
    const int latticeSize = 10;
    const int number_of_sweeps = 100;
    constexpr int steps_per_sweep = latticeSize*latticeSize;
    float temperature = 2;
    float inv_temperature = 1/temperature;
    lattice isingMagnet;

    initializeSimulation(latticeSize, isingMagnet);
    for (int n=0; n < number_of_sweeps; n++){
        for (int m=0; m < steps_per_sweep; m++){
            int rand_x = floor(latticeSize*dist(rd));
            int rand_y = floor(latticeSize*dist(rd));

            int spinFlipCandidate = spinFlipEnergyChange(isingMagnet, rand_x, rand_y);

            if(exp(-inv_temperature*spinFlipCandidate > dist(rd))){
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

//void initializeRender(){}

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