//
// Created by bruno on 01/11/2020.
//

#ifndef CORRIDORALLOCATION_GACORRIDOR_H
#define CORRIDORALLOCATION_GACORRIDOR_H


#include "utils.h"
#include "constants.h"

typedef struct GACorridor {
    double sizeTop;
    double sizeBottom;
    double commCost;
    int indexCut;
    int numRooms;
    int currentNumRooms;
    int currentNumRoomsTop;
    int currentNumRoomsBottom;
    int* corridor;
    double* posMid;

    GACorridor(){}

    explicit GACorridor(int numRooms)  {
        this->commCost = INF;
        this->sizeTop = 0;
        this->sizeBottom = 0;
        this->indexCut = 0;
        this->currentNumRooms = 0;
        this->currentNumRoomsTop = 0;
        this->currentNumRoomsBottom = 0;
        this->numRooms = numRooms;
        this->corridor = new int [numRooms];
        this->posMid = new double [numRooms];
        for(int i = 0; i < numRooms; i ++){
            this->corridor[i] = -1;
            this->posMid[i] = 0;
        }

    }
} GACorridor;


int findSmallestSide(GACorridor * c);
/** roomA, roomB: index in corridor array */
double heuristicaSalas(GACorridor* c, int roomA, int roomB, double* roomsSize, double** commMatrix);
double heuristicaComunicacaoPorTamanho(GACorridor* c, int roomA, int roomB, double* roomsSize, double** commMatrix);

double communicationCost(GACorridor* c, int roomA, int roomB, double* roomsSize, double** commMatrix);
double objectiveFunction(GACorridor* c, double** commMatrix);
void fixRoomPositions(GACorridor* c, double* roomsSize);
void applyMovementSwapRooms(GACorridor* c, int index_a, int index_b, double* roomsSize);

// Given a solution, apply movements to explore neighborhood
void localSearchSwap(GACorridor* c, double* roomsSize, double** commMatrix );

int chooseFirstRoom(int numSalas, double* roomsSize, double** commMatrix, int algorithm);
double getHeuristicCost(GACorridor* c, int roomIndex, double* roomsSize, double** commMatrix);
void moveCorridor(GACorridor* c, int index, bool moveCut);
void insertRoomTop(GACorridor* c, int indexRoom);

void insertRoomBottom(GACorridor* c, int indexRoom);
void removeRoomTop(GACorridor* c);

void removeRoomBottom(GACorridor* c);

void ordenaSalas(GACorridor* c, int side, double* roomsSize, double** commMatrix, bool* isRoomInSolution, std::vector<std::pair <double, int>>* v);

GACorridor buildSolutionRandomized(double* roomsSize, int numRooms, double** commMatrix, double alfa);

GACorridor buildSolutionRandomizedReactive(double* salas, int numSalas, double** commMatrix, int param, double iniAlfa, double fimAlfa, int blocoIt);

void testKnownSolution(double* salas, int numSalas, double** commMatrix);

void buildRandomSolution(GACorridor* c, double* rooms, int numRooms, double** commMatrix);

void initializePopulation(GACorridor* population, double* rooms, int numRooms, double** commMatrix);
int selectElite(GACorridor* currPop, GACorridor* nextPop, int numRooms);
int tournament(GACorridor* currPop);
void crossOver(GACorridor* indivA, GACorridor* indivB, GACorridor* newIndivA, GACorridor* newIndivB, double* rooms);
void mutation(GACorridor* indiv, double* rooms, double** commMatrix);

void printCorridor(GACorridor* c);

void deepCopy(GACorridor *s2, GACorridor *s1);
void setIndexCut(GACorridor* c, double* rooms);

#endif //CORRIDORALLOCATION_GACORRIDOR_H

/*
for(int x = 0; x < 1; x++){
std::cout << "( " << x << " ) " << std::endl;
srand(x + 55);

GACorridor solution(numRooms);

//Corridor solution;

gettimeofday(&tv, 0);
timeStart = (double)tv.tv_sec + 1.0e-6*(double)tv.tv_usec;

int algorithm_select = 3;
switch(algorithm_select){
case GREEDY:
solution = buildSolutionRandomized(rooms, numRooms, commMatrix, 0);
break;
case RANDOMIZED:
solution = buildSolutionRandomized(rooms, numRooms, commMatrix, 0.3);
break;
case REACTIVE:
solution = buildSolutionRandomizedReactive(rooms, numRooms, commMatrix, 11, 0, 1, 5);
break;
default:
std::cout << "Invalid selected algorithm" << std::endl;
return 0;
}


std::cout << objectiveFunction(&solution, commMatrix) << std::endl;

//Corridor local_search_sol = localSearchSwap(rooms, numRooms, commMatrix, &solution);
localSearchSwap(&solution, rooms, commMatrix);

//Corridor local_search_sol = solution;


gettimeofday(&tv, 0);
timeEnd = (double)tv.tv_sec + 1.0e-6*(double)tv.tv_usec;


std::cout << std::endl;
for(int i = 0; i < local_search_sol.roomsTop.size(); i++){
    std::cout << local_search_sol.roomsTop[i].index << " ";
}
std::cout << std::endl;
for(int i = 0; i < local_search_sol.roomsBottom.size(); i++) {
    std::cout << local_search_sol.roomsBottom[i].index << " ";
}
std::cout << std::endl;
std::cout << local_search_sol.commCost << std::endl;
std::cout << timeEnd-timeStart << std::endl;
if(local_search_sol.commCost < best){
    best = local_search_sol.commCost;
}
mean_best += local_search_sol.commCost;
mean_time += timeEnd-timeStart;

fprintf(f_res, "%f,%f\n", local_search_sol.commCost, timeEnd-timeStart);


std::cout << std::endl;
for(int i = 0; i < solution.indexCut; i++){
std::cout << solution.corridor[i] << " ";
}
std::cout << std::endl;
for(int i = solution.indexCut; i < numRooms; i++) {
std::cout << solution.corridor[i] << " ";
}
std::cout << std::endl;
std::cout << solution.commCost << std::endl;
std::cout << timeEnd-timeStart << std::endl;
if(solution.commCost < best){
best = solution.commCost;
}
mean_best += solution.commCost;
mean_time += timeEnd-timeStart;

fprintf(f_res, "%f,%f\n", solution.commCost, timeEnd-timeStart);
}

fprintf(f_res_treated, "%f,%f,%f\n", best, mean_best/30, mean_time/30);


return 0;*/