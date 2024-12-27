#include <iostream>

#ifndef TOPOJUDGE_HPP_
#define TOPOJUDGE_HPP_

using namespace std;


bool IfContainPoint(double node_x, double node_y, double node_w, double obj_x, double obj_y);
bool IfContainMBB(double node_x, double node_y, double node_w, double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax);
bool IfOverlapMBB(double node_x, double node_y, double width, double obj_x_min, double obj_y_min, double obj_x_max, double obj_y_max);
bool IfInterctSeg(double node_x, double node_y, double node_w, double seg_xMin, double seg_yMin, double seg_xMax, double seg_yMax, string seg_type);
bool IfwithinBox(double x, double y, double box_xMin, double box_yMin, double box_xMax, double box_yMax);
bool IfwithinPolygon(double x, double y, double ring_x[], double ring_y[], int len);

#endif