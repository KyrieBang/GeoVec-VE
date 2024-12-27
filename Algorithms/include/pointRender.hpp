#include "indexNode.hpp"

#ifndef _POINTRENDER_HPP_
#define _POINTRENDER_HPP_


string disDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, pointNode *pointRNode);

string dataDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, pointNode *pointRNode, short switch_level);

pointNode * locTileNode(int z, int x, int y, pointNode *seek);

double DiPGPlot(int z, int x, int y, float width, pointNode *RNode);
double DaPGPlot(int z, int x, int y, float width, pointNode *RNode, short transi_level);

#endif