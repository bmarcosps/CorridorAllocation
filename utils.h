//
// Created by bruno on 01/11/2020.
//

#ifndef CORRIDORALLOCATION_UTILS_H
#define CORRIDORALLOCATION_UTILS_H
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

void readInstance(double** rooms, int* numRooms, double*** commMatrix, const char* filename);
double protectedDiv(double a, double b);
int randomIndex(int start, int end, double alfa);

#endif //CORRIDORALLOCATION_UTILS_H
