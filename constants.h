//
// Created by bruno on 01/11/2020.
//

#ifndef CORRIDORALLOCATION_CONSTANTS_H
#define CORRIDORALLOCATION_CONSTANTS_H

#define TOP 1
#define BOTTOM -1

#define GREEDY 1
#define RANDOMIZED 2
#define REACTIVE 3

#define LOWEST_COST 1
#define LOWEST_COST_SIZE 2
#define RANDOM 3

#define INF 9999999999

/** Genetic parameters */
#ifndef NUM_INDIV
#define NUM_INDIV 1000
#endif // NUM_INDIV

#ifndef PROB_CROSS
#define PROB_CROSS 0.7
#endif // PROB_CROSS

#ifndef PROB_MUT
#define PROB_MUT 0.1
#endif // PROB_MUT

#define ELITISM 0.05
#define TOURNAMENT_SIZE 2

#define NUM_GENERATIONS 10000
#define NUM_GENERATIONS_NO_IMPROVE 1000

#define NUM_RUNS 5

#define NUM_EVALUATIONS 2.40e+007

#endif //CORRIDORALLOCATION_CONSTANTS_H
