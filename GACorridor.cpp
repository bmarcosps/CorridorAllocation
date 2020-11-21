//
// Created by bruno on 01/11/2020.
//

#include "GACorridor.h"
/*
std::vector<GARoom> setupRoomList(int numRooms, double* roomsSize){
    std::vector<GARoom> roomList;
    for(int i = 0; i < numRooms; i++){
        roomList.emplace_back(GARoom(roomsSize[i]));
    }
    return roomList;
}
*/

int findSmallestSide(GACorridor * c){
    return (c->sizeBottom > c->sizeTop) ? TOP : BOTTOM;
}
/** roomA, roomB: index in corridor array */
double heuristicaSalas(GACorridor* c, int roomA, int roomB, double* roomsSize, double** commMatrix){
    //double dist = fabs(salaA->posMid - salaB->posMid) * commMatrix[salaA->index][salaB->index];
    // relação comunicação / inverso da distancia
    if(roomsSize[c->corridor[roomA]] < 0) {
        std::cout << "ERRO TAMANHO SALA" << std::endl;
        exit(1);
    }
    return (protectedDiv((roomsSize[c->corridor[roomB]]),fabs(c->posMid[roomA] - c->posMid[roomB]) * commMatrix[c->corridor[roomA]][c->corridor[roomB]]));
}

double heuristicaComunicacaoPorTamanho(GACorridor* c, int roomA, int roomB, double* roomsSize, double** commMatrix){
    double dist = fabs(c->posMid[roomA] - c->posMid[roomB]) * commMatrix[c->corridor[roomA]][c->corridor[roomB]];
    // relação comunicação / tamanho da sala
    return (roomsSize[c->corridor[roomB]])/dist;
}

double communicationCost(GACorridor* c, int roomA, int roomB, double* roomsSize, double** commMatrix){
    double dist = fabs(c->posMid[roomA] - c->posMid[roomB]) * commMatrix[c->corridor[roomA]][c->corridor[roomB]];
    //salaA->cost += dist;
    //salaB->cost += dist;
    return dist;
}

double objectiveFunction(GACorridor* c, double** commMatrix){
    double custo = 0;
    for(int i = 0; i < c->numRooms; i++) {
        for(int j = i+1; j < c->numRooms; j++) {
            custo += fabs(c->posMid[i] - c->posMid[j]) * commMatrix[c->corridor[i]][c->corridor[j]];
        }
    }
    c->commCost = custo;
    return custo;
}

void fixRoomPositions(GACorridor* c, double* roomsSize){
    double pos = 0;
    /*for(int i = 0; i < c->currentNumRooms; i++){
        c->posMid[i] = pos + (roomsSize[c->corridor[i]] / 2);
        pos += roomsSize[c->corridor[i]];

        if(i == c->indexCut-1) {
            pos = 0;
        }
    }*/

    for(int i = 0; i < c->currentNumRoomsTop; i++){
        c->posMid[i] = pos + (roomsSize[c->corridor[i]] / 2);
        pos += roomsSize[c->corridor[i]];
    }
    c->sizeTop = pos;
    pos = 0;
    for(int i = c->numRooms-1; i > c->numRooms - c->currentNumRoomsBottom - 1; i--){
        c->posMid[i] = pos + (roomsSize[c->corridor[i]] / 2);
        pos += roomsSize[c->corridor[i]];
    }
    c->sizeBottom = pos;
}

void applyMovementSwapRooms(GACorridor* c, int index_a, int index_b, double* roomsSize){
    int aux = c->corridor[index_a];
    c->corridor[index_a] = c->corridor[index_b];
    c->corridor[index_b] = aux;
    fixRoomPositions(c, roomsSize);
}

void undoMovementSwapRooms(GACorridor* c, int index_a, int index_b, double* roomsSize){
    int aux = c->corridor[index_a];
    c->corridor[index_a] = c->corridor[index_b];
    c->corridor[index_b] = aux;
    fixRoomPositions(c, roomsSize);
}

void deepCopy(GACorridor *s2, GACorridor *s1)
{
    s2->commCost = s1->commCost;
    s2->currentNumRooms = s1->currentNumRooms;
    s2->currentNumRoomsTop = s1->currentNumRoomsTop;
    s2->currentNumRoomsBottom = s1->currentNumRoomsBottom;
    s2->numRooms = s1->numRooms;
    s2->indexCut = s1->indexCut;
    s2->sizeBottom = s1->sizeBottom;
    s2->sizeTop = s1->sizeTop;
    for(int i = 0; i < s1->numRooms; i++){
        s2->corridor[i] = s1->corridor[i];
        s2->posMid[i] = s1->posMid[i];
    }
}

// Given a solution, apply movements to explore neighborhood
void localSearchSwap(GACorridor* c, double* roomsSize, double** commMatrix ){
    //GACorridor* bestCorridor = new GACorridor(c->numRooms);
    //GACorridor* localBestCorridor = new GACorridor(c->numRooms);

    //(*bestCorridor) = (*c);
    //deepCopy((bestCorridor) , (c));
    int cont = 0;
    int index_a = -1;
    int last_index = -1;
    int best_index_a, best_index_b;
    double best_cost = c->commCost;

    bool improve = false;
    do {
        improve = false;
        //deepCopy((c) , (bestCorridor));
        best_index_a = 0;
        best_index_b = 0;
        for(int i = 0; i < c->numRooms; i++){
            index_a = i;
            for(int j = i; j < c->numRooms; j++){
                if(j != index_a){
                    //cont++;
                    // (*localBestCorridor) = (*c);
                    //deepCopy((localBestCorridor) , (c));
                    applyMovementSwapRooms(c, index_a, j, roomsSize );
                    //setIndexCut(c, roomsSize);
                    objectiveFunction(c, commMatrix);
                    //if(localBestCorridor->commCost < bestCorridor->commCost){
                    if(c->commCost < best_cost){
                        best_cost = c->commCost;
                        best_index_a = i;
                        best_index_b = j;
                        //deepCopy((bestCorridor) , (localBestCorridor));
                        //(*bestCorridor) = (*localBestCorridor);
                        improve = true;
                    }
                    undoMovementSwapRooms(c, index_a, j, roomsSize );
                }
            }
        }
        applyMovementSwapRooms(c, best_index_a, best_index_b, roomsSize );
        c->commCost = best_cost;
        //objectiveFunction(c, commMatrix);

    } while(improve);

    //return (*bestCorridor);

}

int chooseFirstRoom(int numSalas, double* roomsSize, double** commMatrix, int algorithm){
    int index = -1;
    double menorComm = INF;
    double soma = 0;

    switch(algorithm){
        case LOWEST_COST:
            for(int i = 0; i < numSalas; i++) {
                for (int j = 0; j < numSalas; j++) {
                    soma += commMatrix[i][j];
                }

                if(soma < menorComm){
                    menorComm = soma;
                    index = i;
                }
                soma = 0;
            }
            break;
        case LOWEST_COST_SIZE:
            for(int i = 0; i < numSalas; i++) {
                for (int j = 0; j < numSalas; j++) {
                    soma += commMatrix[i][j];
                }
                soma /= roomsSize[i];
                if(soma < menorComm){
                    menorComm = soma;
                    index = i;
                }
                soma = 0;
            }
            break;
        case RANDOM:
            index = randomIndex(0, numSalas-1, 1);
            break;
        default:
            std::cout<< "Invalid algorithm" << std::endl;
            index = -1;
            break;

    }
    return index;

}

double getHeuristicCost(GACorridor* c, int roomIndex, double* roomsSize, double** commMatrix){
    double custo = 0;
    //calcula valor da heuristica em relação a todas as salas alocadas
    for(int i = 0; i < c->numRooms; i++) {
        if(c->corridor[i] != -1)
            custo += heuristicaSalas(c, i, roomIndex, roomsSize, commMatrix);
    }

    return custo;
}

void moveCorridor(GACorridor* c, int index, bool moveCut){
    for(int i = c->numRooms - 1; i > index; i--){
        c->corridor[i] = c->corridor[i-1];
    }
    if(moveCut){
        c->indexCut ++;
    }
}

void insertRoomTop(GACorridor* c, int indexRoom){
    int i;
    //for(i = c->numRooms - 1; i > c->indexCut; i--){
    //    c->corridor[i] = c->corridor[i-1];
    //}
    //c->indexCut++;

    c->corridor[c->currentNumRoomsTop] = indexRoom;
    c->currentNumRooms++;
    c->currentNumRoomsTop++;
}

void insertRoomBottom(GACorridor* c, int indexRoom){
    //c->corridor[c->currentNumRooms] = indexRoom;
    //c->currentNumRooms++;

    c->corridor[c->numRooms - c->currentNumRoomsBottom - 1] = indexRoom;
    c->currentNumRooms++;
    c->currentNumRoomsBottom++;
}

void removeRoomTop(GACorridor* c){
    /*int i;
    for(i = c->indexCut-1; i < c->currentNumRooms-1; i++){
        c->corridor[i] = c->corridor[i+1];
    }
    c->indexCut--;
    c->currentNumRooms--;*/

    int i;

    /*for(i = c->indexCut-1; i < c->currentNumRooms-1; i++){
        c->corridor[i] = c->corridor[i+1];
    }*/
    //c->indexCut--;
    c->corridor[c->currentNumRoomsTop-1] = -1;
    c->currentNumRooms--;
    c->currentNumRoomsTop--;
}

void removeRoomBottom(GACorridor* c){
    c->corridor[c->numRooms - c->currentNumRoomsBottom] = -1;
    c->currentNumRooms--;
    c->currentNumRoomsBottom--;
}

void ordenaSalas(GACorridor* c, int side, double* roomsSize, double** commMatrix, bool* isRoomInSolution, std::vector<std::pair <double, int>>* v){
    double candidateCost = 0;
    for(int i = 0; i < c->numRooms; i++){
        if(!isRoomInSolution[i]) {

            if(side == TOP){
                insertRoomTop(c, i);
            } else {
                insertRoomBottom(c, i);
            }
            fixRoomPositions(c, roomsSize);

            /*
            if(side == TOP){
                candidateCost = getHeuristicCost(c, c->indexCut-1, roomsSize, commMatrix);
            } else {
                candidateCost = getHeuristicCost(c, c->currentNumRooms-1, roomsSize, commMatrix);
            }
            */
            if(side == TOP){
                candidateCost = getHeuristicCost(c, c->currentNumRoomsTop, roomsSize, commMatrix);
            } else {
                candidateCost = getHeuristicCost(c, c->numRooms - c->currentNumRoomsBottom, roomsSize, commMatrix);
            }

            // salva o custo da heuristica do candidato e o indice da sala (no vetor de posicoes)
            (*v).emplace_back(std::make_pair(candidateCost, i));

            if(side == TOP){
                removeRoomTop(c);
            } else {
                removeRoomBottom(c);
            }
        }
    }
    //ordena as roomsSize de acordo com o custo heurístico
    std::sort((*v).begin(), (*v).end());
}

GACorridor buildSolutionRandomized(double* roomsSize, int numRooms, double** commMatrix, double alfa){

    GACorridor *c;
    c = new GACorridor(numRooms);

    // qual lado esta colocando a sala
    int currentSide = TOP;

    // verificação se todas as roomsSize estao alocadas
    int numRoomsInSolution = 0;
    bool* isRoomInSolution = new bool[numRooms];
    for(int i = 0; i < numRooms; i++){
        isRoomInSolution[i] = false;
    }

    // salar ordenadas de acordo com o custo dado pela heurística
    std::vector<std::pair <double, int>> sortedRooms;
    sortedRooms.clear();

    // seleciona primeira sala a entrar na solução
    int index = chooseFirstRoom(numRooms, roomsSize, commMatrix, RANDOM);


    c->corridor[0] = index;
    c->posMid[0] = roomsSize[index] / 2;
    c->sizeTop = roomsSize[index];
    //c->indexCut = 1;
    c->currentNumRooms++;
    c->currentNumRoomsTop++;

    isRoomInSolution[index] = true;
    numRoomsInSolution++;

    while (numRoomsInSolution < numRooms){

        // aloca sempre no lado com menor tamaho total
        currentSide = findSmallestSide(c);

        // calcula o custo heuristico de inserir cada uma das roomsSize ainda nao alocadas
        ordenaSalas(c, currentSide, roomsSize, commMatrix, isRoomInSolution, &sortedRooms);


        //Room* auxRoom = new Room;
        //auxRoom->side = currentSide;

        // escolhe a sala a ser inserida de acordo com alfa
        int newRoomIndex;
        newRoomIndex = randomIndex(0, sortedRooms.size()-1, alfa);
        //std::cout << indiceSala << " ";


        //insere a sala no corredor
        //auxRoom->index = sortedRooms[indiceSala].second;
        newRoomIndex = sortedRooms[newRoomIndex].second;

        if(currentSide == TOP){
            insertRoomTop(c, newRoomIndex);
            c->posMid[c->currentNumRoomsTop-1] = c->sizeTop + roomsSize[newRoomIndex]/2;
            c->sizeTop +=  roomsSize[newRoomIndex];
        } else {
            insertRoomBottom(c, newRoomIndex);
            c->posMid[c->numRooms - c->currentNumRoomsTop] = c->sizeBottom + roomsSize[newRoomIndex]/2;
            c->sizeBottom += roomsSize[newRoomIndex];
        }


        isRoomInSolution[newRoomIndex] = true;
        numRoomsInSolution++;

        sortedRooms.clear();
    }
    c->indexCut = c->currentNumRoomsTop;
/*
    std::cout << std::endl;
    for(int i = 0; i < c->roomsTop.size(); i++){
        std::cout << c->roomsTop[i].index << " ";
    }
    std::cout << std::endl;
    for(int i = 0; i < c->roomsBottom.size(); i++) {
        std::cout << c->roomsBottom[i].index << " ";
    }
    std::cout << std::endl;
*/
    objectiveFunction(c, commMatrix);

    //std::cout << c->commCost << std::endl;

    return (*c);
}

GACorridor buildSolutionRandomized(GACorridor* c, double* roomsSize, int numRooms, double** commMatrix, double alfa){
    c->commCost = INF;
    c->sizeTop = 0;
    c->sizeBottom = 0;
    c->indexCut = 0;
    c->currentNumRooms = 0;
    c->currentNumRoomsTop = 0;
    c->currentNumRoomsBottom = 0;
    //c = new GACorridor(numRooms);

    // qual lado esta colocando a sala
    int currentSide = TOP;

    // verificação se todas as roomsSize estao alocadas
    int numRoomsInSolution = 0;
    bool* isRoomInSolution = new bool[numRooms];
    for(int i = 0; i < numRooms; i++){
        isRoomInSolution[i] = false;
    }

    // salar ordenadas de acordo com o custo dado pela heurística
    std::vector<std::pair <double, int>> sortedRooms;
    sortedRooms.clear();

    // seleciona primeira sala a entrar na solução
    int index = chooseFirstRoom(numRooms, roomsSize, commMatrix, RANDOM);


    c->corridor[0] = index;
    c->posMid[0] = roomsSize[index] / 2;
    c->sizeTop = roomsSize[index];
    //c->indexCut = 1;
    c->currentNumRooms++;
    c->currentNumRoomsTop++;

    isRoomInSolution[index] = true;
    numRoomsInSolution++;

    while (numRoomsInSolution < numRooms){

        // aloca sempre no lado com menor tamaho total
        currentSide = findSmallestSide(c);

        // calcula o custo heuristico de inserir cada uma das roomsSize ainda nao alocadas
        ordenaSalas(c, currentSide, roomsSize, commMatrix, isRoomInSolution, &sortedRooms);


        //Room* auxRoom = new Room;
        //auxRoom->side = currentSide;

        // escolhe a sala a ser inserida de acordo com alfa
        int newRoomIndex;
        newRoomIndex = randomIndex(0, sortedRooms.size()-1, alfa);
        //std::cout << indiceSala << " ";


        //insere a sala no corredor
        //auxRoom->index = sortedRooms[indiceSala].second;
        newRoomIndex = sortedRooms[newRoomIndex].second;

        if(currentSide == TOP){
            insertRoomTop(c, newRoomIndex);
            c->posMid[c->currentNumRoomsTop-1] = c->sizeTop + roomsSize[newRoomIndex]/2;
            c->sizeTop +=  roomsSize[newRoomIndex];
        } else {
            insertRoomBottom(c, newRoomIndex);
            c->posMid[c->numRooms - c->currentNumRoomsTop] = c->sizeBottom + roomsSize[newRoomIndex]/2;
            c->sizeBottom += roomsSize[newRoomIndex];
        }


        isRoomInSolution[newRoomIndex] = true;
        numRoomsInSolution++;

        sortedRooms.clear();
    }
    c->indexCut = c->currentNumRoomsTop;

    objectiveFunction(c, commMatrix);

    return (*c);
}

GACorridor buildSolutionRandomizedReactive(double* salas, int numSalas, double** commMatrix, int param, double iniAlfa, double fimAlfa, int blocoIt){

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


    GACorridor *c = new GACorridor(numSalas);
    GACorridor *bestCorridor;
    bestCorridor = new GACorridor(numSalas);
    bestCorridor->commCost = INF;

    for(int j = 0; j < param; j++){
        *c = buildSolutionRandomized(salas, numSalas, commMatrix, alfas[j]);
        somaSolucoes[j] = somaSolucoes[j] + c->commCost;

        if(c->commCost < bestCorridor->commCost){
            (*bestCorridor) = (*c);
        }

    }

    for(int iteracao = 1; iteracao < 50; iteracao++){
        //imprimeInformacoes(alfas, probAlfa, usoAlfa, param);
        std::cout << iteracao << " ";

        float prob = ((float)rand()/(float)(RAND_MAX+1));
        int indiceAtual = 0;
        while (prob - probAlfa[indiceAtual] > 0){
            prob = prob - probAlfa[indiceAtual];
            indiceAtual++;
        }
        usoAlfa[indiceAtual]++;

        *c = buildSolutionRandomized( salas, numSalas, commMatrix, alfas[indiceAtual]);
        std::cout << (*c).commCost << std::endl;

        localSearchSwap(c, salas, commMatrix);
        //*c = localSearch(salas, numSalas, commMatrix, c);

        somaSolucoes[indiceAtual] += somaSolucoes[indiceAtual] + c->commCost;

        if(c->commCost < bestCorridor->commCost){
            (*bestCorridor) = (*c);
        }
        std::cout << (*bestCorridor).commCost << std::endl;
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

void testKnownSolution(double* salas, int numSalas, double** commMatrix){
    GACorridor* c = new GACorridor(numSalas);
    c->corridor[0]=5;
    c->corridor[1]=1;
    c->corridor[2]=9;
    c->corridor[3]=0;
    c->corridor[4]=2;

    c->corridor[9]=7;
    c->corridor[8]=3;
    c->corridor[7]=4;
    c->corridor[6]=6;
    c->corridor[5]=8;

    //c->indexCut = 5;

    c->currentNumRooms = 10;
    c->currentNumRoomsTop = 5;
    c->currentNumRoomsBottom = 5;

    fixRoomPositions(c, salas);
    for(int i = 0; i < c->numRooms; i++){
        std::cout << salas[c->corridor[i]] << " ";
    }
    std::cout << std::endl;
    for(int i = 0; i < c->numRooms; i++){
        std::cout << c->posMid[i] << " ";
    }
    objectiveFunction(c, commMatrix);
    std::cout << c->commCost << std::endl;
}


void printCorridor(GACorridor* c){
    std::cout << std::endl;
    for(int i = 0; i < c->currentNumRoomsTop; i++){
        std::cout << c->corridor[i] << " ";
    }
    std::cout << std::endl;
    for(int i = c->numRooms -1; i > c->numRooms-c->currentNumRoomsBottom - 1; i--) {
        std::cout << c->corridor[i] << " ";
    }
    std::cout << std::endl;
    std::cout << c->commCost << std::endl;
}

void setIndexCut(GACorridor* c, double* rooms){
    int i = 0;
    int j = c->numRooms-1;
    c->currentNumRoomsTop = 1;
    c->currentNumRoomsBottom = 1;

    double topSize = rooms[c->corridor[i++]];
    double bottomSize = rooms[c->corridor[j--]];

    while (i <= j){
        if(topSize > bottomSize){
            c->currentNumRoomsBottom++;
            bottomSize += rooms[c->corridor[j--]];
        } else {
            c->currentNumRoomsTop++;
            topSize += rooms[c->corridor[i++]];
        }
    }
    c->sizeTop = topSize;
    c->sizeBottom = bottomSize;

    /*float probMiddle = ((float)rand()/(float)(RAND_MAX+1));
    if(probMiddle <= 0.8) {
        c->sizeTop = topSize;
        c->sizeBottom = bottomSize;
    } else {
        probMiddle = ((float)rand()/(float)(RAND_MAX+1));
        if(probMiddle <= 0.5) {
            c->sizeTop = topSize + 1;
            c->sizeBottom = bottomSize - 1;
        } else {
            c->sizeTop = topSize - 1;
            c->sizeBottom = bottomSize + 1;
        }
    }
*/
}

void setOptimalIndexCut(GACorridor* c, double* rooms, double** commMatrix){
    int bestTopSize = 0;
    double bestObj = INF;
    int i = 0;
    int j = c->numRooms-1;
    for(int k = 1; k < c->numRooms-1; k++){
        c->currentNumRoomsTop = k;
        c->currentNumRoomsBottom = c->numRooms - k;
        fixRoomPositions(c, rooms);
        objectiveFunction(c, commMatrix);
        if(c->commCost < bestObj) {
            bestObj = c->commCost;
            bestTopSize = k;
        }
    }
    c->currentNumRoomsTop = bestTopSize;
    c->currentNumRoomsBottom = c->numRooms - bestTopSize;
    fixRoomPositions(c, rooms);

    //c->currentNumRoomsTop = 1;
    //c->currentNumRoomsBottom = c->numRooms - 1;

    //double topSize = rooms[c->corridor[i++]];
    //double bottomSize = rooms[c->corridor[j--]];
/*
    while (i <= j){
        if(topSize > bottomSize){
            c->currentNumRoomsBottom++;
            bottomSize += rooms[c->corridor[j--]];
        } else {
            c->currentNumRoomsTop++;
            topSize += rooms[c->corridor[i++]];
        }
    }
*/
}

void buildRandomSolution(GACorridor* c, double* rooms, int numRooms, double** commMatrix){
    //std::cout << "Building random Solution" << std::endl;
    std::vector<int> v(numRooms);
    std::iota(v.begin(), v.end(), 0);
    for(int i = 0; i < numRooms; i++){
        int index = randomIndex(0, v.size()-1, 1);
        c->corridor[i] = v[index];
        v.erase(v.begin()+index);
    }
    c->currentNumRooms = numRooms;
    setIndexCut(c, rooms);
    fixRoomPositions(c, rooms);
    /*
    int x = c->sizeTop;
    int y = c->sizeBottom;
    setOptimalIndexCut(c, rooms, commMatrix);
    if(c->sizeTop != x){
        std::cout << x << " " << y << std::endl;
        std::cout << c->sizeTop << " " << c->sizeBottom << std::endl;

    }
*/


}


void initializePopulation(GACorridor* population, double* rooms, int numRooms, double** commMatrix){
    //5% constructive
    int i = 0;
    //for(; i < 0.05 * NUM_INDIV; i++){
        //population[i] = buildSolutionRandomizedReactive(rooms, numRooms, commMatrix, 11, 0, 1, 5);
        //buildSolutionRandomized(&population[i] , rooms, numRooms, commMatrix, 0.5);
        //localSearchSwap(&population[i], rooms, commMatrix);
    //}
    for(i = 0; i < NUM_INDIV; i++){
        buildRandomSolution(&population[i], rooms, numRooms, commMatrix);
    }
    //95% random
}

void sortElite(GACorridor* pop, int k, int numRooms){
    GACorridor* aux = new GACorridor(numRooms);
    int i, j;
    int bestIndex;

    for(i = 0; i < k; i++){
        bestIndex = i;
        for(j = i+1; j < NUM_INDIV; j++){
            if(pop[j].commCost < pop[bestIndex].commCost){
                bestIndex = j;
            }
        }
        //aux = pop[i];

        deepCopy(aux, &pop[i]);
        deepCopy(&pop[i], &pop[bestIndex]);
        deepCopy(&pop[bestIndex], aux);


        //pop[i] = pop[bestIndex];
        //pop[bestIndex] = aux;
    }
    delete aux;
}


int selectElite(GACorridor* currPop, GACorridor* nextPop, int numRooms){
    int i;
    int cont = (int)(ELITISM * NUM_INDIV);//(ELITISMO*NUM_INDIV);

    if(cont % 2 != 0)
        cont++;

    sortElite(currPop, cont, numRooms);
    for(i = 0; i < cont; i++){
        deepCopy(&nextPop[i], &currPop[i]);
        //popFutura[i] = popAtual[i];
    }

    return cont;
}

int tournament(GACorridor* currPop){
    int bestIndex = randomIndex(0, NUM_INDIV-1, 1);
    int index;
    for(int i = 0; i < TOURNAMENT_SIZE-1; i++){
        index = randomIndex(0, NUM_INDIV-1, 1);
        if(currPop[index].commCost < currPop[bestIndex].commCost){
            bestIndex = index;
        }
    }
    return bestIndex;
}

void crossOverOX(GACorridor* indivA, GACorridor* indivB, GACorridor* newIndiv){
    int i = randomIndex(0, indivA->numRooms-1, 1)% indivA->numRooms;
    int j = randomIndex(0, indivA->numRooms-1, 1) % indivA->numRooms;
    while( abs(i-j) > (int)(0.75 * indivA->numRooms) ){
        i = randomIndex(0, indivA->numRooms-1, 1)% indivA->numRooms;
        j = randomIndex(0, indivA->numRooms-1, 1) % indivA->numRooms;
    }
    if(j < i) {
        int aux = i;
        i = j;
        j = aux;
    }
    newIndiv->currentNumRooms = 0;

    bool* isRoomInSolution = new bool[indivA->numRooms];
    for(int x = 0; x < indivA->numRooms; x++){
        isRoomInSolution[x] = false;
    }

    int k;
    for(k = i; k <=j; k++){
        newIndiv->corridor[k] = indivA->corridor[k];
        isRoomInSolution[indivA->corridor[k]] = true;
        newIndiv->currentNumRooms++;
    }
    int l = k;
    while(newIndiv->currentNumRooms != indivA->numRooms){
        if( l == indivA->numRooms) l = 0;
        if( k == indivA->numRooms) k = 0;
        if(!isRoomInSolution[indivB->corridor[l]]){
            newIndiv->corridor[k++] = indivB->corridor[l++];
            newIndiv->currentNumRooms++;
        } else {
            l++;
        }
    }

}

void crossOver(GACorridor* indivA, GACorridor* indivB, GACorridor* newIndivA, GACorridor* newIndivB, double* rooms){
    crossOverOX(indivA, indivB, newIndivA);
    setIndexCut(newIndivA, rooms);
    fixRoomPositions(newIndivA, rooms);

    crossOverOX(indivB, indivA, newIndivB);
    setIndexCut(newIndivB, rooms);
    fixRoomPositions(newIndivB, rooms);
}



void mutation(GACorridor* c, double* rooms, double** commMatrix){
   // std::cout << "before " << indiv->commCost << std::endl;


    float mutProb = ((float)rand()/(float)(RAND_MAX+1));
    int i = randomIndex(0, c->numRooms, 1)% c->numRooms;
    int j = randomIndex(0, c->numRooms-1, 1) % c->numRooms;
    while(i == j)
        j = randomIndex(0, c->numRooms-1, 1) % c->numRooms;

    applyMovementSwapRooms(c, i, j, rooms );
    setIndexCut(c, rooms);
    fixRoomPositions(c, rooms);

    if (mutProb > 0.90) {
        localSearchSwap(c, rooms, commMatrix);
        //objectiveFunction(c, commMatrix);
    }
    //std::cout << "after " << indiv->commCost << std::endl;
}