#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <math.h>
#include <algorithm>

#define TOP 1
#define BOTTOM -1

#define ALGORITHM 2

#define GREEDY 1
#define RANDOMIZED 2

#define INF 999999999

typedef struct Room {
    unsigned int index;
    double posStart, posMid, posEnd;

    Room(int index, double x, double y, double z){
        this->index = index;
        this->posStart = x;
        this->posMid = y;
        this->posEnd = z;
    }

    Room(){
        this->index = -1;
        this->posStart = 0;
        this->posMid = 0;
        this->posEnd = 0;
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

void readInstance(double** rooms, int* numRooms, double*** commMatrix, char* filename){
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

double heuristicaSalas(Room* salaA, Room* salaB, double** commMatrix){
    double dist = fabs(salaA->posMid - salaB->posMid) * commMatrix[salaA->index][salaB->index];
    return dist + (protectedDiv(1,fabs(salaA->posMid - salaB->posMid)) * commMatrix[salaA->index][salaB->index]);
}

double distanciaSalas(Room* salaA, Room* salaB, double** commMatrix){
    double dist = fabs(salaA->posMid - salaB->posMid) * commMatrix[salaA->index][salaB->index];
    return dist;
}

double custoCommTotal(Corridor* c, double** commMatrix){
    double custo = 0;
    for(int i = 0; i < c->roomsTop.size(); i++){
        for(int j = i+1; j < c->roomsTop.size(); j++){
            custo += distanciaSalas(&c->roomsTop[i], &c->roomsTop[j], commMatrix);
        }

        for(int k = 0; k < c->roomsBottom.size(); k++){
            custo += distanciaSalas(&c->roomsTop[i], &c->roomsBottom[k], commMatrix);
        }
    }

    for(int i = 0; i < c->roomsBottom.size(); i++) {
        for (int j = i + 1; j < c->roomsBottom.size(); j++) {
            custo += distanciaSalas(&c->roomsBottom[i], &c->roomsBottom[j], commMatrix);
        }
    }
    return custo;
}

int primeiraSala(int numSalas, double** commMatrix){
    double menorComm = INF;
    int menorIndex = -1;
    double soma = 0;
    for(int i = 0; i < numSalas; i++) {
        for (int j = 0; j < numSalas; j++) {
            soma += commMatrix[i][j];
        }
        if(soma < menorComm){
            menorComm = soma;
            menorIndex = i;
        }
        soma = 0;
    }
    return menorIndex;
}


int primeiraSala2(int numSalas, double** commMatrix){
    double maiorComm = -1;
    int maiorIndex = -1;
    double soma = 0;
    for(int i = 0; i < numSalas; i++) {
        for (int j = 0; j < numSalas; j++) {
            soma += commMatrix[i][j];
        }
        if(soma > maiorComm){
            maiorComm = soma;
            maiorIndex = i;
        }
        soma = 0;
    }
    return maiorIndex;
}

double calculaDistancia(Room* s, Corridor* c, double** commMatrix){
    double custo = 0;
    for(auto & i : c->roomsTop){
        custo += heuristicaSalas(&i, s, commMatrix);
    }

    for(int i = 0; i < c->roomsBottom.size(); i++) {
        custo += heuristicaSalas(&c->roomsBottom[i], s, commMatrix);
    }
    return custo;
}

void ordenaSalas(int corredor, Corridor* c, double* salas, double** commMatrix, int numSalas, bool* salaSolucao, std::vector<std::pair <double, int>>* v){
    for(int i = 0; i < numSalas; i++){
        if(!salaSolucao[i]) {
            Room s;
            s.index = i;

            if(corredor == TOP){
                if(c->roomsTop.empty()){
                    s.posStart = 0;
                } else {
                    s.posStart = c->roomsTop.back().posEnd;
                }
            } else {
                if(c->roomsBottom.empty()){
                    s.posStart = 0;
                } else {
                    s.posStart = c->roomsBottom.back().posEnd;
                }
            }

            s.posEnd = s.posStart + salas[i];
            s.posMid = s.posStart + (salas[i] / 2);

            (*v).emplace_back(std::make_pair(calculaDistancia(&s,c,commMatrix), i));
        }
    }
    std::sort((*v).begin(), (*v).end());
}


int randomIndex(int start, int end, double alfa){
    return start + (rand() % ((int)(alfa * end)+1));
}
int findSmallestSide(Corridor* c){
    return (c->sizeBottom > c->sizeTop) ? TOP : BOTTOM;
}

Corridor buildSolution(double* salas, int numSalas, double** commMatrix, double alfa){

    Corridor *c;
    c = new Corridor();

    int currentSide = TOP;

    int numRoomsInSolution = 0;
    bool* isRoomInSolution = new bool[numSalas];

    std::vector<std::pair <double, int>> sortedRooms;
    sortedRooms.clear();

    for(int i = 0; i < numSalas; i++){
        isRoomInSolution[i] = false;
    }

    int index = primeiraSala(numSalas, commMatrix);


    Room* s = new Room();

    s->index = index;
    s->posStart = 0;
    s->posEnd = salas[index];
    s->posMid = salas[index] / 2;
    c->roomsTop.push_back(*s);
    c->sizeBottom = salas[index];

    isRoomInSolution[index] = true;
    numRoomsInSolution++;
    while (numRoomsInSolution < numSalas){
        //currentSide *= -1;
        currentSide = findSmallestSide(c);
        ordenaSalas(currentSide, c, salas, commMatrix, numSalas, isRoomInSolution, &sortedRooms);

        Room* s = new Room;
        unsigned int indiceSala;

        switch(ALGORITHM){
            case GREEDY:
                indiceSala = 0;
                break;
            case RANDOMIZED:
                indiceSala = randomIndex(0, sortedRooms.size(), alfa);
                break;
            default:
                break;
        }


        s->index = sortedRooms[indiceSala].second;
        if(currentSide == TOP){
            if(c->roomsTop.empty()){
                s->posStart = 0;
            } else {
                s->posStart = c->roomsTop.back().posEnd;
            }
        } else {
            if(c->roomsBottom.empty()){
                s->posStart = 0;
            } else {
                s->posStart = c->roomsBottom.back().posEnd;
            }
        }
        s->posEnd = s->posStart + salas[s->index];
        s->posMid = s->posStart + (salas[s->index] / 2);

        if(currentSide == TOP){
            c->roomsTop.push_back(*s);
        } else {
            c->roomsBottom.push_back(*s);
        }

        //std::cout << "-" << s->posStart << "-" <<std::endl;
        //std::cout <<  sortedRooms[0].second <<std::endl;
        //std::cout <<  sortedRooms[0].first <<std::endl;

        c->commCost += sortedRooms[indiceSala].first;
        isRoomInSolution[s->index] = true;
        numRoomsInSolution++;
        sortedRooms.clear();
    }
    for(int i = 0; i < c->roomsTop.size(); i++){
        std::cout << c->roomsTop[i].index << " ";
    }
    std::cout << std::endl;
    for(int i = 0; i < c->roomsBottom.size(); i++) {
        std::cout << c->roomsBottom[i].index << " ";
    }
    std::cout << std::endl;

    std::cout << c->commCost;

    return (*c);
}



int main(int argc, char** argv) {
    //8
    srand(8);

    char* datasetFile = argv[1];

    int numRooms;
    double* rooms;
    double** commMatrix;

    readInstance(&rooms, &numRooms, &commMatrix, datasetFile);

    for(int i = 0; i < (numRooms); i++){
        for(int j = 0; j<(numRooms); j++){
            std::cout << commMatrix[i][j]<< " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    Corridor c = buildSolution(rooms, numRooms, commMatrix, 0.3);

    std::cout << std::endl;

    std::cout << custoCommTotal(&c, commMatrix) << std::endl;
    return 0;
}
