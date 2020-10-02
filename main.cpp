#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <math.h>
#include <algorithm>

#define CIMA 1
#define BAIXO -1

#define ALGORITMO 2
#define GULOSO 1
#define RANDOMIZADO 2

#define INF 999999999

typedef struct Sala{
    unsigned int indice;
    double posIni, posMeio, posFim;

    Sala(int index, double x, double y, double z){
        indice = index;
        posIni = x;
        posMeio = y;
        posFim = z;
    }

    Sala(){
        indice = -1;
        posIni = 0;
        posMeio = 0;
        posFim = 0;
    }
} Sala;

typedef struct Corredor{
    std::vector<Sala> filaCima;
    std::vector<Sala> filaBaixo;
    double custoComunicacao;
    Corredor(){
        custoComunicacao = 0;
    }
} Corredor;

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

double heuristicaSalas(Sala* salaA, Sala* salaB, double** commMatrix){
    double dist = fabs(salaA->posMeio - salaB->posMeio) * commMatrix[salaA->indice][salaB->indice];
    return dist + (protectedDiv(1,fabs(salaA->posMeio - salaB->posMeio)) * commMatrix[salaA->indice][salaB->indice]);
}

double distanciaSalas(Sala* salaA, Sala* salaB, double** commMatrix){
    double dist = fabs(salaA->posMeio - salaB->posMeio) * commMatrix[salaA->indice][salaB->indice];
    return dist;
}

double custoCommTotal(Corredor* c, double** commMatrix){
    double custo = 0;
    for(int i = 0; i < c->filaCima.size(); i++){
        for(int j = i+1; j < c->filaCima.size(); j++){
            custo += distanciaSalas(&c->filaCima[i], &c->filaCima[j], commMatrix);
        }

        for(int k = 0; k < c->filaBaixo.size(); k++){
            custo += distanciaSalas(&c->filaCima[i], &c->filaBaixo[k], commMatrix);
        }
    }

    for(int i = 0; i < c->filaBaixo.size(); i++) {
        for (int j = i + 1; j < c->filaBaixo.size(); j++) {
            custo += distanciaSalas(&c->filaBaixo[i], &c->filaBaixo[j], commMatrix);
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

double calculaDistancia(Sala* s, Corredor* c, double** commMatrix){
    double custo = 0;
    for(int i = 0; i < c->filaCima.size(); i++){
        custo += heuristicaSalas(&c->filaCima[i], s, commMatrix);
    }

    for(int i = 0; i < c->filaBaixo.size(); i++) {
        custo += heuristicaSalas(&c->filaBaixo[i], s, commMatrix);
    }
    return custo;
}

void ordenaSalas(int corredor, Corredor* c, double* salas, double** commMatrix, int numSalas, bool* salaSolucao, std::vector<std::pair <double, int>>* v){
    for(int i = 0; i < numSalas; i++){
        if(!salaSolucao[i]) {
            Sala s;
            s.indice = i;

            if(corredor == CIMA){
                if(c->filaCima.empty()){
                    s.posIni = 0;
                } else {
                    s.posIni = c->filaCima.back().posFim;
                }
            } else {
                if(c->filaBaixo.empty()){
                    s.posIni = 0;
                } else {
                    s.posIni = c->filaBaixo.back().posFim;
                }
            }

            s.posFim = s.posIni + salas[i];
            s.posMeio = s.posIni + (salas[i]/2);

            (*v).emplace_back(std::make_pair(calculaDistancia(&s,c,commMatrix), i));
        }
    }
    std::sort((*v).begin(), (*v).end());
}


int randomIndex(int start, int end, double alfa){
    return start + (rand() % ((int)(alfa * end)+1));
}

Corredor construirSolucao(double* salas, int numSalas, double** commMatrix, double alfa){

    Corredor *c;
    c = new Corredor();

    int lado = CIMA;

    int numSolucao = 0;
    bool* salaSolucao = new bool[numSalas];

    std::vector<std::pair <double, int>> salasOrdenadas;
    salasOrdenadas.clear();

    for(int i = 0; i < numSalas; i++){
        salaSolucao[i] = false;
    }

    int index = primeiraSala(numSalas, commMatrix);


    Sala* s = new Sala();
    s->indice = index;
    s->posIni = 0;
    s->posFim = salas[index];
    s->posMeio = salas[index]/2;
    c->filaCima.push_back(*s);

    salaSolucao[index] = true;
    numSolucao++;
    while (numSolucao < numSalas){
        lado *= -1;
        ordenaSalas(lado, c, salas, commMatrix, numSalas, salaSolucao, &salasOrdenadas);

        Sala* s = new Sala;
        unsigned int indiceSala;

        switch(ALGORITMO){
            case GULOSO:
                indiceSala = 0;
                break;
            case RANDOMIZADO:
                indiceSala = randomIndex(0, salasOrdenadas.size(), alfa);
                break;
            default:
                break;
        }


        s->indice = salasOrdenadas[indiceSala].second;
        if(lado == CIMA){
            if(c->filaCima.empty()){
                s->posIni = 0;
            } else {
                s->posIni = c->filaCima.back().posFim;
            }
        } else {
            if(c->filaBaixo.empty()){
                s->posIni = 0;
            } else {
                s->posIni = c->filaBaixo.back().posFim;
            }
        }
        s->posFim = s->posIni + salas[s->indice];
        s->posMeio = s->posIni + (salas[s->indice]/2);

        if(lado == CIMA){
            c->filaCima.push_back(*s);
        } else {
            c->filaBaixo.push_back(*s);
        }

        //std::cout << "-" << s->posIni << "-" <<std::endl;
        //std::cout <<  salasOrdenadas[0].second <<std::endl;
        //std::cout <<  salasOrdenadas[0].first <<std::endl;

        c->custoComunicacao += salasOrdenadas[indiceSala].first;
        salaSolucao[s->indice] = true;
        numSolucao++;
        salasOrdenadas.clear();
    }
    for(int i = 0; i < c->filaCima.size(); i++){
        std::cout << c->filaCima[i].indice << " ";
    }
    std::cout << std::endl;
    for(int i = 0; i < c->filaBaixo.size(); i++) {
        std::cout << c->filaBaixo[i].indice << " ";
    }
    std::cout << std::endl;

    std::cout << c->custoComunicacao;

    return (*c);
}



int main(int argc, char** argv) {
    //8
    srand(8);
    char* datasetFile = argv[1];
    double* salas;
    int numSalas;
    double** commMatrix;

    readInstance(&salas, &numSalas, &commMatrix, datasetFile);

    for(int i = 0; i < (numSalas); i++){
        for(int j = 0; j<(numSalas); j++){
            std::cout << commMatrix[i][j]<< " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    Corredor c = construirSolucao( salas, numSalas, commMatrix, 0.3);

    std::cout << std::endl;

    std::cout << custoCommTotal(&c, commMatrix) << std::endl;
    return 0;
}
