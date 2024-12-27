#include "indexNode.hpp"

#ifndef _LINESTRINGRENDER_HPP_
#define _LINESTRINGRENDER_HPP_


string disDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, linestringNode *linestringRNode);

string dataDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, linestringNode *linestringRNode, short switch_level);

linestringNode * locTileNode(int z, int x, int y, linestringNode *seek);

double DiPGPlot(int z, int x, int y, float width, linestringNode *RNode);
double DaPGPlot(int z, int x, int y, float width, linestringNode *RNode, short transi_level);

#endif