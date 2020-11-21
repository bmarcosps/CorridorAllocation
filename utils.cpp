//
// Created by bruno on 01/11/2020.
//

#include "utils.h"
#include "GACorridor.h"

void readInstance(double** rooms, int* numRooms, double*** commMatrix, const char* filename){
    std::fstream arq;
    std::string line;

    printf("Lendo Dados do Arquivo... %s\n",filename);
    arq.open(filename, std::fstream::in);

    arq >> (*numRooms);

    std::cout << (*numRooms) << " ";

    (*rooms) = new double[*numRooms];
    (*commMatrix) = new double* [*numRooms];
    for(int i = 0; i < (*numRooms); i++){
        (*commMatrix)[i] = new double[*numRooms];
    }

    arq >> line;
    std::istringstream iss(line);
    std::string token;
    for(int i = 0; i<(*numRooms); i++){
        std::getline(iss, token, ',');
        (*rooms)[i] = std::stof(token);
        std::cout << (*rooms)[i] << " ";
    }

    std::cout << std::endl;
    for(int i = 0; i < (*numRooms); i++){
        arq >> line;
        std::istringstream iss2(line);
        for(int j = 0; j<(*numRooms); j++){
            std::getline(iss2, token, ',');
            (*commMatrix)[i][j] = std::stof(token);
            //std::cout << commMatrix[i][j]<< " ";
        }
        //std::cout << std::endl;
    }
}

double protectedDiv(double a, double b){
    if (b == 0){
        return 1;
    } else {
        return a/b;
    }
}

int randomIndex(int start, int end, double alfa){
    return start + (rand() % ((int)(alfa * end) + 1));
}

