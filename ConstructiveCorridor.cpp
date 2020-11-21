//
// Created by bruno on 01/11/2020.
//

#include "ConstructiveCorridor.h"


int findSmallestSide(Corridor* c){
    return (c->sizeBottom > c->sizeTop) ? TOP : BOTTOM;
}

double heuristicaSalas(Room* salaA, Room* salaB, double** commMatrix){
    //double dist = fabs(salaA->posMid - salaB->posMid) * commMatrix[salaA->index][salaB->index];
    // relação comunicação / inverso da distancia
    if(salaB->posEnd - salaB->posStart < 0) {
        std::cout << "ERRO TAMANHO SALA" << std::endl;
        exit(1);
    }
    return (protectedDiv((salaB->posEnd - salaB->posStart),fabs(salaA->posMid - salaB->posMid) * commMatrix[salaA->index][salaB->index]));
}

double heuristicaComunicacaoPorTamanho(Room* salaA, Room* salaB, double** commMatrix){
    double dist = fabs(salaA->posMid - salaB->posMid) * commMatrix[salaA->index][salaB->index];
    // relação comunicação / tamanho da sala
    return (salaB->posEnd - salaB->posStart)/dist;
}

double communicationCost(Room* salaA, Room* salaB, double** commMatrix){
    double dist = fabs(salaA->posMid - salaB->posMid) * commMatrix[salaA->index][salaB->index];
    salaA->cost += dist;
    salaB->cost += dist;
    return dist;
}

double objectiveFunction(Corridor* c, double** commMatrix){
    double custo = 0;
    for(int i = 0; i < c->roomsTop.size(); i++){
        for(int j = i+1; j < c->roomsTop.size(); j++){
            custo += communicationCost(&c->roomsTop[i], &c->roomsTop[j], commMatrix);
        }

        for(int k = 0; k < c->roomsBottom.size(); k++){
            custo += communicationCost(&c->roomsTop[i], &c->roomsBottom[k], commMatrix);
        }
    }

    for(int i = 0; i < c->roomsBottom.size(); i++) {
        for (int j = i + 1; j < c->roomsBottom.size(); j++) {
            custo += communicationCost(&c->roomsBottom[i], &c->roomsBottom[j], commMatrix);
        }
    }
    c->commCost = custo;
    return custo;
}

Room findRoom(Corridor* c, int roomIndex){
    double custo = 0;
    for(int i = 0; i < c->roomsTop.size(); i++) {
        if(c->roomsTop[i].index == roomIndex){
            return c->roomsTop[i];
        }
    }
    for(int i = 0; i < c->roomsBottom.size(); i++) {
        if(c->roomsBottom[i].index == roomIndex){
            return c->roomsBottom[i];
        }
    }


    //return NULL;
}

void fixRoomPositions(Corridor* c, double* salas){
    double pos = 0;
    for(int i = 0; i < c->roomsTop.size(); i++) {
        c->roomsTop.at(i).posStart = pos;

        c->roomsTop.at(i).posEnd = pos + salas[c->roomsTop.at(i).index];
        c->roomsTop.at(i).posMid = pos + (salas[c->roomsTop.at(i).index] / 2);
        c->roomsTop.at(i).indexAtCorridor = i;
        c->roomsTop.at(i).cost = 0;
        c->roomsTop.at(i).side = TOP;
        pos += salas[c->roomsTop.at(i).index];
    }

    pos = 0;
    for(int i = 0; i < c->roomsBottom.size(); i++) {
        c->roomsBottom.at(i).posStart = pos;

        c->roomsBottom.at(i).posEnd = pos + salas[c->roomsBottom.at(i).index];
        c->roomsBottom.at(i).posMid = pos + (salas[c->roomsBottom.at(i).index] / 2);
        c->roomsBottom.at(i).indexAtCorridor = i;
        c->roomsBottom.at(i).cost = 0;
        c->roomsBottom.at(i).side = BOTTOM;
        pos += salas[c->roomsBottom.at(i).index];
    }
}

void applyMovementMoveRoom(Corridor* c, double* salas, int numSalas, int index_a, int index_b){
    Room r_a = findRoom(c, index_a);
    Room r_b;
    if(index_b < numSalas) {
        r_b = findRoom(c, index_b);
        if(r_a.side == r_b.side){
            if(r_a.side == TOP){
                c->roomsTop.erase(c->roomsTop.begin() + r_a.indexAtCorridor); // at(r_a.indexAtCorridor) = r_b;
                c->roomsTop.insert(c->roomsTop.begin() + r_b.indexAtCorridor, r_a);
            } else {
                c->roomsBottom.erase(c->roomsBottom.begin() + r_a.indexAtCorridor); // at(r_a.indexAtCorridor) = r_b;
                c->roomsBottom.insert(c->roomsBottom.begin() + r_b.indexAtCorridor, r_a);
            }
        } else {
            if(r_a.side == TOP){
                c->roomsTop.erase(c->roomsTop.begin() + r_a.indexAtCorridor); // at(r_a.indexAtCorridor) = r_b;
                c->roomsBottom.insert(c->roomsBottom.begin() + r_b.indexAtCorridor, r_a);
            } else {
                c->roomsBottom.erase(c->roomsBottom.begin() + r_a.indexAtCorridor); // at(r_a.indexAtCorridor) = r_b;
                c->roomsTop.insert(c->roomsTop.begin() + r_b.indexAtCorridor, r_a);
            }
        }
    } else {
        if(r_a.side == TOP){
            c->roomsTop.erase(c->roomsTop.begin() + r_a.indexAtCorridor); // at(r_a.indexAtCorridor) = r_b;
            c->roomsBottom.insert(c->roomsBottom.end() , r_a);
        } else {
            c->roomsBottom.erase(c->roomsBottom.begin() + r_a.indexAtCorridor); // at(r_a.indexAtCorridor) = r_b;
            c->roomsTop.insert(c->roomsTop.end(), r_a);
        }
    }


    fixRoomPositions(c, salas);
}

void applyMovementSwapRooms(Corridor* c, double* salas, int numSalas, int index_a, int index_b){
    Room r_a = findRoom(c, index_a);
    Room r_b = findRoom(c, index_b);

    if(r_a.side == r_b.side){
        if(r_a.side == TOP){

            c->roomsTop.at(r_a.indexAtCorridor) = r_b;
            c->roomsTop.at(r_b.indexAtCorridor) = r_a;
        } else {

            c->roomsBottom.at(r_a.indexAtCorridor) = r_b;
            c->roomsBottom.at(r_b.indexAtCorridor) = r_a;
        }
    } else {
        if(r_a.side == TOP){

            c->roomsTop.at(r_a.indexAtCorridor) = r_b;
            c->roomsBottom.at(r_b.indexAtCorridor) = r_a;
        } else {

            c->roomsBottom.at(r_a.indexAtCorridor) = r_b;
            c->roomsTop.at(r_b.indexAtCorridor) = r_a;
        }
    }

    fixRoomPositions(c, salas);
}

int getWorstCost(Corridor* c, double alfa){
    double cost = 0;
    int index = 0;
    std::vector<std::pair <double, int>> v;
    for(int i = 0; i < c->roomsTop.size(); i++) {
        (v).emplace_back(std::make_pair(c->roomsTop[i].cost, c->roomsTop[i].index));
    }

    for(int i = 0; i < c->roomsBottom.size(); i++) {
        (v).emplace_back(std::make_pair(c->roomsBottom[i].cost, c->roomsBottom[i].index));
    }

    std::sort((v).begin(), (v).end(), std::greater <>());

    int indiceSala = randomIndex(0, v.size()-1, alfa);

    return v[indiceSala].second;
}

// Given a solution, apply movements to explore neighborhood
Corridor localSearch(double* salas, int numSalas, double** commMatrix, Corridor* c){
    int iterations = 1000;
    //int i = 0;
    Corridor* bestCorridor = new Corridor();
    Corridor* localBestCorridor = new Corridor();

    (*bestCorridor) = (*c);

    int index_a = -1;
    int last_index = -1;
    do {
        (*c) = (*bestCorridor);
        //for(int i = 0; i < numSalas; i++){
        index_a = getWorstCost(c, 0.4);


        for(int j = 0; j <= numSalas; j++){
            (*localBestCorridor) = (*c);

            if(j != index_a){
                applyMovementMoveRoom(localBestCorridor, salas, numSalas, index_a, j);
                objectiveFunction(localBestCorridor, commMatrix);
                if(localBestCorridor->commCost < bestCorridor->commCost){
                    std::cout << "Improvement" << std::endl;
                    (*bestCorridor) = (*localBestCorridor);
                }
            }
        }
        //}

        iterations--;
    } while(bestCorridor->commCost - c->commCost != 0 && iterations > 0);

    return (*bestCorridor);

}

Corridor localSearchSwap(double* salas, int numSalas, double** commMatrix, Corridor* c){
    Corridor* bestCorridor = new Corridor();
    Corridor* localBestCorridor = new Corridor();
    struct timeval tv;
    double timeStart, timeEnd;
    (*bestCorridor) = (*c);
    int cont = 0;
    int index_a = -1;
    int last_index = -1;
    bool improve = false;
    do {
        improve = false;
        (*c) = (*bestCorridor);
        //std::vector<int> indexes(numSalas);
        //std::iota(indexes.begin(), indexes.end(), 0);
        //gettimeofday(&tv, 0);
        //timeStart = (double)tv.tv_sec + 1.0e-6*(double)tv.tv_usec;
        for(int i = 0; i < numSalas; i++){

            index_a = i;
            for(int j = i; j < numSalas; j++){
                if(j != index_a){
                    //cont++;
                    (*localBestCorridor) = (*c);
                    applyMovementSwapRooms(localBestCorridor, salas, numSalas, index_a, j);
                    objectiveFunction(localBestCorridor, commMatrix);

                    if(localBestCorridor->commCost < bestCorridor->commCost){
                        (*bestCorridor) = (*localBestCorridor);
                        improve = true;
                    }
                }
            }
        }
        //gettimeofday(&tv, 0);
        //timeEnd = (double)tv.tv_sec + 1.0e-6*(double)tv.tv_usec;
        //std::cout << timeEnd-timeStart << " " << cont << std::endl;

    } while(improve);

    return (*bestCorridor);

}

void addRoomToCorridor(double* salas, Corridor* c, Room* r, int side){
    if(side == TOP){
        if(c->roomsTop.empty()){
            r->posStart = 0;
        } else {
            r->posStart = c->roomsTop.back().posEnd;
        }
    } else {
        if(c->roomsBottom.empty()){
            r->posStart = 0;
        } else {
            r->posStart = c->roomsBottom.back().posEnd;
        }
    }
    r->posEnd = r->posStart + salas[r->index];
    r->posMid = r->posStart + (salas[r->index] / 2);

    if(side == TOP){
        c->roomsTop.push_back(*r);
        c->sizeTop += salas[r->index];
    } else {
        c->roomsBottom.push_back(*r);
        c->sizeBottom += salas[r->index];
    }

}

int chooseFirstRoom(int numSalas, double* salas, double** commMatrix, int algorithm){
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
                soma /= salas[i];
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

double getHeuristicCost(Room* s, Corridor* c, double** commMatrix){
    double custo = 0;
    //calcula valor da heuristica em relação a todas as salas alocadas
    for(int i = 0; i < c->roomsTop.size(); i++) {
        //custo += heuristicaComunicacaoPorTamanho(&c->roomsTop[i], s, commMatrix);
        custo += heuristicaSalas(&c->roomsTop[i], s, commMatrix);
    }

    for(int i = 0; i < c->roomsBottom.size(); i++) {
        //custo += heuristicaComunicacaoPorTamanho(&c->roomsBottom[i], s, commMatrix);
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

            double custoCandidato = getHeuristicCost(&s, c, commMatrix);

            // salva o custo da heuristica do candidato e o indice da sala (no vetor de posicoes)
            (*v).emplace_back(std::make_pair(custoCandidato, i));
        }
    }

    //ordena as salas de acordo com o custo heurístico
    std::sort((*v).begin(), (*v).end());
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
    int index = chooseFirstRoom(numSalas, salas, commMatrix, RANDOM);


    Room* s = new Room();

    s->index = index;
    s->posStart = 0;
    s->posEnd = salas[index];
    s->posMid = salas[index] / 2;
    s->indexAtCorridor = 0;
    s->side = TOP;
    c->roomsTop.push_back(*s);
    c->sizeTop = salas[index];

    isRoomInSolution[index] = true;
    numRoomsInSolution++;
    while (numRoomsInSolution < numSalas){

        // aloca sempre no lado com menor tamaho total
        currentSide = findSmallestSide(c);

        // calcula o custo heuristico de inserir cada uma das salas ainda nao alocadas
        ordenaSalas(currentSide, c, salas, commMatrix, numSalas, isRoomInSolution, &sortedRooms);


        Room* auxRoom = new Room;
        auxRoom->side = currentSide;

        // escolhe a sala a ser inserida de acordo com alfa
        unsigned int indiceSala;
        indiceSala = randomIndex(0, sortedRooms.size()-1, alfa);
        //std::cout << indiceSala << " ";

        //insere a sala no corredor
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
            auxRoom->indexAtCorridor = c->roomsTop.size();

            c->roomsTop.push_back(*auxRoom);
            c->sizeTop += salas[auxRoom->index];
        } else {
            auxRoom->indexAtCorridor = c->roomsBottom.size();
            c->roomsBottom.push_back(*auxRoom);
            c->sizeBottom += salas[auxRoom->index];
        }


        isRoomInSolution[auxRoom->index] = true;
        numRoomsInSolution++;
        sortedRooms.clear();
    }
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
        *c = localSearchSwap(salas, numSalas, commMatrix, c);
        *c = localSearch(salas, numSalas, commMatrix, c);

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
    Corridor* c = new Corridor();


    Room* r = new Room();
    r->index = 5;
    addRoomToCorridor(salas, c, r, TOP);
    delete r;
    r = new Room();
    r->index = 1;
    addRoomToCorridor(salas, c, r, TOP);
    delete r;
    r = new Room();
    r->index = 9;
    addRoomToCorridor(salas, c, r, TOP);
    delete r;
    r = new Room();
    r->index = 0;
    addRoomToCorridor(salas, c, r, TOP);
    delete r;
    r = new Room();
    r->index = 2;
    addRoomToCorridor(salas, c, r, TOP);
    delete r;
    r = new Room();
    r->index = 7;
    addRoomToCorridor(salas, c, r, BOTTOM);
    delete r;
    r = new Room();
    r->index = 3;
    addRoomToCorridor(salas, c, r, BOTTOM);
    delete r;
    r = new Room();
    r->index = 4;
    addRoomToCorridor(salas, c, r, BOTTOM);
    delete r;
    r = new Room();
    r->index = 6;
    addRoomToCorridor(salas, c, r, BOTTOM);
    delete r;
    r = new Room();
    r->index = 8;
    addRoomToCorridor(salas, c, r, BOTTOM);
    delete r;

    objectiveFunction(c, commMatrix);
    std::cout << c->commCost << std::endl;
}