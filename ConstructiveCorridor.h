//
// Created by bruno on 01/11/2020.
//

#ifndef CORRIDORALLOCATION_CONSTRUCTIVECORRIDOR_H
#define CORRIDORALLOCATION_CONSTRUCTIVECORRIDOR_H

#include "utils.h"
#include "constants.h"

typedef struct Room {
    //índice original da sala no vetor de tamanhos
    unsigned int index;
    unsigned int indexAtCorridor;
    int side;
    double posStart, posMid, posEnd;
    double cost;

    Room(int index, double x, double y, double z, int indexAtCorridor, int side, double cost){
        this->index = index;
        this->posStart = x;
        this->posMid = y;
        this->posEnd = z;
        this->indexAtCorridor = indexAtCorridor;
        this->side = side;
        this->cost = cost;
    }

    Room(){
        this->index = -1;
        this->posStart = 0;
        this->posMid = 0;
        this->posEnd = 0;
        this->indexAtCorridor = -1;
        this->side = -1;
        this->cost = 0;
    }
} Room;

typedef struct Corridor{
    std::vector<Room> roomsTop;
    std::vector<Room> roomsBottom;
    double commCost;
    double sizeTop, sizeBottom;

    Corridor(){
        commCost = 0;
        sizeTop = 0;
        sizeBottom = 0;
    }
} Corridor;

int findSmallestSide(Corridor* c);

double heuristicaSalas(Room* salaA, Room* salaB, double** commMatrix);

double heuristicaComunicacaoPorTamanho(Room* salaA, Room* salaB, double** commMatrix);

double communicationCost(Room* salaA, Room* salaB, double** commMatrix);

double objectiveFunction(Corridor* c, double** commMatrix);

Room findRoom(Corridor* c, int roomIndex);

void fixRoomPositions(Corridor* c, double* salas);

void applyMovementMoveRoom(Corridor* c, double* salas, int numSalas, int index_a, int index_b);

void applyMovementSwapRooms(Corridor* c, double* salas, int numSalas, int index_a, int index_b);

int getWorstCost(Corridor* c, double alfa);

// Given a solution, apply movements to explore neighborhood
Corridor localSearch(double* salas, int numSalas, double** commMatrix, Corridor* c);

Corridor localSearchSwap(double* salas, int numSalas, double** commMatrix, Corridor* c);

void addRoomToCorridor(double* salas, Corridor* c, Room* r, int side);

int chooseFirstRoom(int numSalas, double* roomsSize, double** commMatrix, int algorithm);

double getHeuristicCost(Room* s, Corridor* c, double** commMatrix);

void ordenaSalas(int corredor, Corridor* c, double* salas, double** commMatrix, int numSalas, bool* salaSolucao, std::vector<std::pair <double, int>>* v);

Corridor buildSolutionRandomized(double* roomsSize, int numRooms, double** commMatrix, double alfa);

Corridor buildSolutionRandomizedReactive(double* salas, int numSalas, double** commMatrix, int param, double iniAlfa, double fimAlfa, int blocoIt);

void testKnownSolution(double* salas, int numSalas, double** commMatrix);

#endif //CORRIDORALLOCATION_CONSTRUCTIVECORRIDOR_H
