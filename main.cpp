#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <numeric>

#include "utils.h"
#include "constants.h"
#include "GACorridor.h"


int main(int argc, char** argv) {
    //8
    //char* datasetFile = argv[1];
    struct timeval tv;
    std::string intanceFolder = ".\\instances\\";
    std::string datasetFile = argv[1];

    std::string filename = ".\\results\\2result_"+ datasetFile +".txt";
    std::string filename_res = ".\\results\\1result_"+ datasetFile +".txt";

    FILE *f_res = fopen(filename.c_str(), "w");
    FILE *f_res_treated = fopen(filename_res.c_str(), "w");


    double mean_best = 0;
    double mean_time = 0;
    double timeStart, timeEnd;

    int numRooms;
    double* rooms;
    double** commMatrix;

    readInstance(&rooms, &numRooms, &commMatrix, intanceFolder.append(datasetFile).c_str());

    for(int i = 0; i < (numRooms); i++){
        for(int j = 0; j<(numRooms); j++){
            std::cout << commMatrix[i][j]<< " ";
        }
        std::cout << std::endl;
    }
    //testKnownSolution(rooms, numRooms, commMatrix);
    std::cout << std::endl;

    GACorridor best_corridor(numRooms);
    GACorridor overall_best(numRooms);

    GACorridor *currentPopulation;
    currentPopulation = new GACorridor[NUM_INDIV];

    GACorridor *nextPopulation;
    nextPopulation = new GACorridor[NUM_INDIV];
    for(int i = 0; i < NUM_INDIV; i++){
        currentPopulation[i] = GACorridor(numRooms);
        nextPopulation[i] = GACorridor(numRooms);
    }



    for(int x = 0; x < NUM_RUNS; x++){
        std::cout << "( " << x << " ) " << std::endl;
        srand(x + 33);

        gettimeofday(&tv, 0);
        timeStart = (double)tv.tv_sec + 1.0e-6*(double)tv.tv_usec;

        //initialize population
        initializePopulation(currentPopulation, rooms, numRooms, commMatrix);

        //evaluate population
        for(int i = 0; i < NUM_INDIV; i++){
            objectiveFunction(&currentPopulation[i], commMatrix);
        }
        int indexA;
        int indexB;
        bool improved;
        int cont_gent;
        float prob_mut = PROB_MUT;

        improved = false;
        double best = INF;
        cont_gent = 0;

        for( int k = 0; k < NUM_GENERATIONS; k++){
            std::cout << "Generation " << k  << std::endl;


            int newIndividuals = selectElite(currentPopulation, nextPopulation, numRooms);
            if(currentPopulation[0].commCost < best){
                cont_gent = 0;
                prob_mut = PROB_MUT;
                deepCopy(&best_corridor, &currentPopulation[0]);
                best = best_corridor.commCost;
            } else if(cont_gent == NUM_GENERATIONS_NO_IMPROVE){
                break;
            } else {

                cont_gent++;
                if(cont_gent % 100 == 0){
                    prob_mut += 0.05;
                }
                std::cout << "( " << cont_gent << " ) " << std::endl;

            }

            std::cout << "Best = " << best << std::endl;
            int j;
            for(j = newIndividuals; j < 0.9 * NUM_INDIV; j +=2){

                //selection
                //int id = (j/2)-(newIndividuals/2);
                indexA = tournament(currentPopulation);
                indexB = tournament(currentPopulation);
                while (indexA == indexB)
                    indexB = tournament(currentPopulation);
                //printCorridor(&currentPopulation[indexA]);
                //printCorridor(&currentPopulation[indexB]);

                //popFutura[j] = popAtual[indice1];
                //popFutura[j+1] = popAtual[indice2];

                float crossProb = ((float)rand()/(float)(RAND_MAX+1));//randomProb();


                if(crossProb <= PROB_CROSS){
                    crossOver(&currentPopulation[indexA], &currentPopulation[indexB], &nextPopulation[j], &nextPopulation[j+1], rooms);
                } else {
                    //nextPopulation[j] = buildRandomSolution(rooms,numRooms);
                    //nextPopulation[j+1] = buildRandomSolution(rooms,numRooms);
                    deepCopy(&nextPopulation[j], &currentPopulation[indexA]);
                    deepCopy(&nextPopulation[j+1], &currentPopulation[indexB]);
                }

                float mutProb = ((float)rand()/(float)(RAND_MAX+1));//randomProb();
                if(mutProb <= prob_mut){
                    //std::cout << "Mutation" << std::endl;
                    mutation(&nextPopulation[j], rooms, commMatrix);
                }
                mutProb = ((float)rand()/(float)(RAND_MAX+1));
                if(mutProb <= prob_mut){
                    //std::cout << "Mutation" << std::endl;
                    mutation(&nextPopulation[j+1], rooms, commMatrix);
                }
                //std::cout << "new " << std::endl;
                //printCorridor(&nextPopulation[j]);
                //printCorridor(&nextPopulation[j+1]);
            }
            for(; j < NUM_INDIV; j ++){
                 buildRandomSolution(&nextPopulation[j], rooms,numRooms,commMatrix);
            }
            for(int i = newIndividuals; i < NUM_INDIV; i++){
                deepCopy(&currentPopulation[i], &nextPopulation[i]);
                objectiveFunction(&currentPopulation[i], commMatrix);
                //std::cout << " " << currentPopulation[i].commCost;
            }
            std::cout << std::endl;

        }
        mean_best += best_corridor.commCost;

        localSearchSwap(&best_corridor, rooms, commMatrix);

        if(best_corridor.commCost < overall_best.commCost){
            deepCopy(&overall_best, &best_corridor);
        }


        gettimeofday(&tv, 0);
        timeEnd = (double)tv.tv_sec + 1.0e-6*(double)tv.tv_usec;
        std::cout << timeEnd-timeStart << std::endl;

        //mean_best += solution.commCost;
        mean_time += timeEnd-timeStart;
        printCorridor(&best_corridor);

        fprintf(f_res, "%f, %f\n", best_corridor.commCost, timeEnd-timeStart);
    }

    printCorridor(&overall_best);

    fprintf(f_res_treated, "%f,%f,%f\n", overall_best.commCost, mean_best/NUM_RUNS, mean_time/NUM_RUNS);


    return 0;
}
