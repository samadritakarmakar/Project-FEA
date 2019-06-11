#ifndef FEMTOOLS_H
#define FEMTOOLS_H
#include "LoadGaussFile.h"

void GetNodePostions(umat &NodePositons, const umat &ElementNodes, int vectorLevel);
void GetNodePostions(umat &NodePositons, const umat &ElementNodes, int vectorLevel, int originalVctrLvl);

#endif
