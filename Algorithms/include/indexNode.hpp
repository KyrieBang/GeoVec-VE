#include <iostream>
#include <fstream>
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>

using namespace std;
namespace bg	= boost::geometry;
namespace bgi	= boost::geometry::index;
namespace bgm	= boost::geometry::model;

#ifndef _INDEXNODE_HPP_
#define _INDEXNODE_HPP_



typedef bgm::d2::point_xy<double> POINT;
typedef bgm::segment<POINT>	SEG;
typedef bgm::box<POINT> Box;
typedef vector<double> XY_VEC;
typedef bgi::rtree<POINT, bgi::quadratic<8>> pointRtree;
typedef pair<Box, XY_VEC> BOX_PAIR;
typedef bgi::rtree<BOX_PAIR, bgi::quadratic<8>> linestringRtree;
typedef bgi::rtree<BOX_PAIR, bgi::quadratic<8>> polygonRtree;



enum nodeType{
    ROOT,
    LU,
    RU,
    LB,
    RB,
};



class pointNode{
    public:
        // node type
        nodeType node_type;
        // node BBox
        double node_x;
        double node_y;
        double node_w;
        // node attributes
        unsigned short node_level;
        unsigned int node_fcount;
        // node ptrs
        pointNode *parent;
        pointNode *LUNode;
        pointNode *RUNode;
        pointNode *LBNode;
        pointNode *RBNode;
        pointRtree Rtree;
        // node functions
        pointNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, pointNode *_parent);
        ~pointNode();
        void insertObj(double obj_x, double obj_y, short transi_level);
        void locRtree(vector<pointRtree> &Rtree_list, short transi_level);
};

class linestringNode{
    public:
        // node type
        nodeType node_type;
        // node BBox
        double node_x; 
        double node_y;
        double node_w;
        // node attributes
        unsigned short node_level;
        unsigned int node_fcount;
        unsigned int node_scount;
        double node_flen;
        // node ptrs
        linestringNode *parent;
        linestringNode *LUNode;
        linestringNode *RUNode;
        linestringNode *LBNode;
        linestringNode *RBNode;
        linestringRtree Rtree;
        // node functions
        linestringNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, linestringNode *_parent);
        ~linestringNode();
        linestringNode *insertTMBB(double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax, int seg_count, double line_len, short transi_level);
        void insertSeg(double seg_xMin, double seg_yMin, double seg_xMax, double seg_yMax, string seg_type, double line_len, short transi_level);
        void insertFeature(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short transi_level, BOX_PAIR line);
        void locRtree(vector<linestringRtree> &Rtree_list, short transi_level);
};

class polygonNode{
    public:
        // node type
        nodeType node_type;
        char topo_type;
        // node BBox
        double node_x; 
        double node_y;
        double node_w;
        // node attributes
        unsigned short node_level;
        unsigned int node_fcount;
        double node_flen;
        double node_farea;
        // node ptrs
        polygonNode *parent;
        polygonNode *LUNode;
        polygonNode *RUNode;
        polygonNode *LBNode;
        polygonNode *RBNode;
        polygonRtree Rtree;
        // node functions
        polygonNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, polygonNode *_parent);
        ~polygonNode();
        polygonNode *insertTMBB(double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax, int edge_count, double obj_len, double obj_area, short transi_level);
        void insertEdge(double edge_xMin, double edge_yMin, double edge_xMax, double edge_yMax, string edge_type, double obj_len, short transi_level);
        void insertFeature(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short transi_level, double ring_x[], double ring_y[], int len, double obj_area, BOX_PAIR poly);
        void locRtree(vector<polygonRtree> &Rtree_list, short transi_level);
        void searchTile(int z, int x, int y, int count, int transi_level, vector<int> &tile_id);
};


#endif