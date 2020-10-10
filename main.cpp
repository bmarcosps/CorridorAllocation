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

#define ALGORITHM 3

#define GREEDY 1
#define RANDOMIZED 2
#define REACTIVE 3


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
    // relação comunicação / inverso da distancia
    return dist + (protectedDiv(1,fabs(salaA->posMid - salaB->posMid)) * commMatrix[salaA->index][salaB->index]);
}

double heuristicaSalas2(Room* salaA, Room* salaB, double** commMatrix){
    double dist = fabs(salaA->posMid - salaB->posMid) * commMatrix[salaA->index][salaB->index];
    // relação comunicação / tamanho da sala
    return dist/(salaB->posEnd - salaB->posStart);
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
    c->commCost = custo;
    return custo;
}

int primeiraSala(int numSalas, double* salas, double** commMatrix, bool menor_custo){
    if(menor_custo){
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
    } else {
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
}


double calculaCustoHeuristica(Room* s, Corridor* c, double** commMatrix){
    double custo = 0;
    for(int i = 0; i < c->roomsTop.size(); i++) {
        custo += heuristicaSalas2(&c->roomsTop[i], s, commMatrix);
    }

    for(int i = 0; i < c->roomsBottom.size(); i++) {
        custo += heuristicaSalas2(&c->roomsBottom[i], s, commMatrix);
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

            double custoCandidato = calculaCustoHeuristica(&s, c, commMatrix);

            (*v).emplace_back(std::make_pair(custoCandidato, i));
        }
    }
    std::sort((*v).begin(), (*v).end());
}


int randomIndex(int start, int end, double alfa){
    return start + (rand() % ((int)(alfa * end) + 1));
}
int findSmallestSide(Corridor* c){
    return (c->sizeBottom > c->sizeTop) ? TOP : BOTTOM;
}

Corridor buildSolutionRandomized(double* salas, int numSalas, double** commMatrix, double alfa){

    Corridor *c;
    c = new Corridor();

    // qual lado esta colocando a sala
    int currentSide = TOP;

    // verificação se todas as salas estao alocadas
    int numRoomsInSolution = 0;
    bool* isRoomInSolution = new bool[numSalas];
    for(int i = 0; i < numSalas; i++){
        isRoomInSolution[i] = false;
    }

    // salar ordenadas de acordo com o custo dado pela heurística
    std::vector<std::pair <double, int>> sortedRooms;
    sortedRooms.clear();

    // seleciona primeira sala a entrar na solução
    // true = primeira é a de menor custo de comunicação total
    // false = primeira é a de maior custo de comunicação total
    int index = primeiraSala(numSalas, salas, commMatrix, true);


    Room* s = new Room();

    s->index = index;
    s->posStart = 0;
    s->posEnd = salas[index];
    s->posMid = salas[index] / 2;
    c->roomsTop.push_back(*s);
    c->sizeTop = salas[index];

    isRoomInSolution[index] = true;
    numRoomsInSolution++;
    while (numRoomsInSolution < numSalas){
        //currentSide *= -1;
        // aloca sempre no lado com menor tamaho total
        currentSide = findSmallestSide(c);
        //std::cout << currentSide << std::endl;
        ordenaSalas(currentSide, c, salas, commMatrix, numSalas, isRoomInSolution, &sortedRooms);

        Room* auxRoom = new Room;
        unsigned int indiceSala;

        indiceSala = randomIndex(0, sortedRooms.size()-1, alfa);

        auxRoom->index = sortedRooms[indiceSala].second;

        if(currentSide == TOP){
            if(c->roomsTop.empty()){
                auxRoom->posStart = 0;
            } else {
                auxRoom->posStart = c->roomsTop.back().posEnd;
            }
        } else {
            if(c->roomsBottom.empty()){
                auxRoom->posStart = 0;
            } else {
                auxRoom->posStart = c->roomsBottom.back().posEnd;
            }
        }
        auxRoom->posEnd = auxRoom->posStart + salas[auxRoom->index];
        auxRoom->posMid = auxRoom->posStart + (salas[auxRoom->index] / 2);

        if(currentSide == TOP){
            c->roomsTop.push_back(*auxRoom);
            c->sizeTop += salas[auxRoom->index];
        } else {
            c->roomsBottom.push_back(*auxRoom);
            c->sizeBottom += salas[auxRoom->index];
        }

        //std::cout << "-" << s->posStart << "-" <<std::endl;
        //std::cout <<  sortedRooms[0].second <<std::endl;
        //std::cout <<  sortedRooms[0].first <<std::endl;

        //c->commCost += sortedRooms[indiceSala].first;
        isRoomInSolution[auxRoom->index] = true;
        numRoomsInSolution++;
        sortedRooms.clear();
    }

    std::cout << std::endl;
    for(int i = 0; i < c->roomsTop.size(); i++){
        std::cout << c->roomsTop[i].index << " ";
    }
    std::cout << std::endl;
    for(int i = 0; i < c->roomsBottom.size(); i++) {
        std::cout << c->roomsBottom[i].index << " ";
    }
    std::cout << std::endl;
    custoCommTotal(c, commMatrix);
    if(c->commCost == 1086){
        std::cout<< "aaa";
    }
    std::cout << c->commCost;

    return (*c);
}

Corridor buildSolutionRandomizedReactive(double* salas, int numSalas, double** commMatrix, int param, double iniAlfa, double fimAlfa, int blocoIt){

    double alfas[param];
    double probAlfa[param];
    double somaSolucoes[param];
    int usoAlfa[param];
    double media[param];
    double Q[param];
    double sumQ = 0;

    for(int u = 0; u < param; u++){
        alfas[u] = iniAlfa + (float)u *((float)(fimAlfa-iniAlfa)/(float)(param-1));
        probAlfa[u] = 1/(float)param;
        usoAlfa[u] = 1;
        somaSolucoes[u] = 0;
        media[u] = 0;
        Q[u] = 0;
    }


    Corridor *c = new Corridor;
    Corridor *bestCorridor;
    bestCorridor = new Corridor();
    bestCorridor->commCost = INF;

    for(int j = 0; j < param; j++){
        *c = buildSolutionRandomized(salas, numSalas, commMatrix, alfas[j]);
        somaSolucoes[j] = somaSolucoes[j] + c->commCost;

        if(c->commCost < bestCorridor->commCost){
            (*bestCorridor) = (*c);
        }

    }

    for(int iteracao = 1; iteracao < 10000; iteracao++){
        //imprimeInformacoes(alfas, probAlfa, usoAlfa, param);


        float prob = ((float)rand()/(float)(RAND_MAX+1));
        int indiceAtual = 0;
        while (prob - probAlfa[indiceAtual] > 0){
            prob = prob - probAlfa[indiceAtual];
            indiceAtual++;
        }
        usoAlfa[indiceAtual]++;
        //saida << "*******************************" << endl;
        //saida << "Iteracao " << iteracao << endl;
        // <<  "Alfa escolhido:  " << alfas[indiceAtual] << endl;
        //cout << "Alfa = " << alfas[indiceAtual] << endl;
        //realiza busca construtiva
        *c = buildSolutionRandomized( salas, numSalas, commMatrix, alfas[indiceAtual]);
        //saida << "Número de vértices: " << r[0] << endl;
        //saida << "Soma dos pesos: " << r[1] << endl;
        //saida << "*******************************" << endl;
        //cout << "r[0] = " << r[0] << endl;
        //cout << "r[1] = " << r[1] << endl;

        somaSolucoes[indiceAtual] += somaSolucoes[indiceAtual] + c->commCost;

        if(c->commCost < bestCorridor->commCost){
            (*bestCorridor) = (*c);
        }
        //delete[] r;
        //free(r);
        //r = NULL;

        if(iteracao%blocoIt == 0){

            //atualiza probabilidades
            for(int u = 0; u < param; u++){

                media[u] = somaSolucoes[u]/usoAlfa[u];
                Q[u] = pow((bestCorridor->commCost/media[u]), 2);
                //cout <<"qiparcial = " << Q[u] << endl;
                sumQ += Q[u];
            }
            //cout <<"SUMQI = " << sumQ << endl;

            for(int j = 0; j < param; j++){
                probAlfa[j] = Q[j] / sumQ;

            }
            sumQ = 0;
        }
    }
    return *bestCorridor;

}


Corridor buildSolution(double* salas, int numSalas, double** commMatrix){

    Corridor *c;
    c = new Corridor();

    // qual lado esta colocando a sala
    int currentSide = TOP;

    // verificação se todas as salas estao alocadas
    int numRoomsInSolution = 0;
    bool* isRoomInSolution = new bool[numSalas];
    for(int i = 0; i < numSalas; i++){
        isRoomInSolution[i] = false;
    }

    // salar ordenadas de acordo com o custo dado pela heurística
    std::vector<std::pair <double, int>> sortedRooms;
    sortedRooms.clear();

    // seleciona primeira sala a entrar na solução
    // true = primeira é a de menor custo de comunicação total
    // false = primeira é a de maior custo de comunicação total
    int index = primeiraSala(numSalas, salas, commMatrix, true);


    Room* s = new Room();

    s->index = index;
    s->posStart = 0;
    s->posEnd = salas[index];
    s->posMid = salas[index] / 2;
    c->roomsTop.push_back(*s);
    c->sizeTop = salas[index];

    isRoomInSolution[index] = true;
    numRoomsInSolution++;
    while (numRoomsInSolution < numSalas){
        //currentSide *= -1;
        // aloca sempre no lado com menor tamaho total
        currentSide = findSmallestSide(c);

        ordenaSalas(currentSide, c, salas, commMatrix, numSalas, isRoomInSolution, &sortedRooms);

        Room* auxRoom = new Room;
        unsigned int indiceSala;

        indiceSala = 0;

        auxRoom->index = sortedRooms[indiceSala].second;
        if(currentSide == TOP){
            if(c->roomsTop.empty()){
                auxRoom->posStart = 0;
            } else {
                auxRoom->posStart = c->roomsTop.back().posEnd;
            }
        } else {
            if(c->roomsBottom.empty()){
                auxRoom->posStart = 0;
            } else {
                auxRoom->posStart = c->roomsBottom.back().posEnd;
            }
        }
        auxRoom->posEnd = auxRoom->posStart + salas[auxRoom->index];
        auxRoom->posMid = auxRoom->posStart + (salas[auxRoom->index] / 2);

        if(currentSide == TOP){
            c->roomsTop.push_back(*auxRoom);
            c->sizeTop = auxRoom->posEnd;
        } else {
            c->roomsBottom.push_back(*auxRoom);
            c->sizeBottom = auxRoom->posEnd;
        }

        //std::cout << "-" << s->posStart << "-" <<std::endl;
        //std::cout <<  sortedRooms[0].second <<std::endl;
        //std::cout <<  sortedRooms[0].first <<std::endl;

        c->commCost += sortedRooms[indiceSala].first;
        isRoomInSolution[auxRoom->index] = true;
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
    srand(10);

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
    Corridor solution;

    switch(ALGORITHM){
        case GREEDY:
            solution = buildSolution(rooms, numRooms, commMatrix);
            break;
        case RANDOMIZED:
            solution = buildSolutionRandomized(rooms, numRooms, commMatrix, 0.3);
            break;
        case REACTIVE:
            solution = buildSolutionRandomizedReactive(rooms, numRooms, commMatrix, 11, 0, 1, 10);
            break;
    }


    std::cout << std::endl;

    std::cout << custoCommTotal(&solution, commMatrix) << std::endl;
    return 0;
}
