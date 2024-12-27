#include "indexNode.hpp"

#ifndef _POLYGONRENDER_HPP_
#define _POLYGONRENDER_HPP_


string disDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, polygonNode *polygonRNode);

string dataDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, polygonNode *polygonRNode, short switch_level);

polygonNode * locTileNode(int z, int x, int y, polygonNode *seek);

double DiPGPlot(int z, int x, int y, float width, polygonNode *RNode);
double DaPGPlot(int z, int x, int y, float width, polygonNode *RNode, short transi_level);


#endif