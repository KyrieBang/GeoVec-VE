#include "indexNode.hpp"
#include "topoJudge.hpp"


// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// ***************************************************    point node functions    *************************************************************** //
// ---------------------------------------------------------------------------------------------------------------------------------------------- //
pointNode::pointNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, pointNode *_parent){
    this->node_type = _node_type;
    this->node_x = _node_x;
    this->node_y = _node_y;
    this->node_w = _node_w;
    this->node_level = _node_level;
    this->node_fcount = 0;
    this->parent = _parent;
    this->LUNode = nullptr;
    this->RUNode = nullptr;
    this->LBNode = nullptr;
    this->RBNode = nullptr;
}


pointNode::~pointNode()
{
    delete LUNode;
    delete RUNode;
    delete LBNode;
    delete RBNode;
    delete parent;
    LUNode = nullptr;
    RUNode = nullptr;
    LBNode = nullptr;
    RBNode = nullptr;
    parent = nullptr;
}


// inert point obj
void pointNode::insertObj(double obj_x, double obj_y, short transi_level){
    if (node_level == transi_level + 1){
        Rtree.insert(POINT(obj_x, obj_y));
    }
    
    if (node_level == transi_level + 8){
        return;
    }

    // LU sub-region
    if (IfContainPoint(node_x, node_y + node_w/2, node_w/2, obj_x, obj_y)){
        if (!LUNode){
            LUNode = new pointNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        LUNode->node_fcount++;
        LUNode->insertObj(obj_x, obj_y, transi_level);
        return;
    }
    // RU sub-region
    else if (IfContainPoint(node_x + node_w/2, node_y + node_w/2, node_w/2, obj_x, obj_y)){
        if (!RUNode){
            RUNode = new pointNode(RU, node_x + node_w/2, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        RUNode->node_fcount++;
        RUNode->insertObj(obj_x, obj_y, transi_level);
        return;
    }
    // LB sub-region
    else if (IfContainPoint(node_x, node_y, node_w/2, obj_x, obj_y)){
        if (!LBNode){
            LBNode = new pointNode(LB, node_x, node_y, node_w/2, node_level + 1, this);
        }
        LBNode->node_fcount++;
        LBNode->insertObj(obj_x, obj_y, transi_level);
        return;
    }
    // RB sub-reion
    else if (IfContainPoint(node_x + node_w/2, node_y, node_w/2, obj_x, obj_y)){
        if (!RBNode){
            RBNode = new pointNode(RB, node_x + node_w/2, node_y, node_w/2, node_level + 1, this);
        }
        RBNode->node_fcount++;
        RBNode->insertObj(obj_x, obj_y, transi_level);
        return;
    }
}


// loc the R-tree list of a node
void pointNode::locRtree(vector<pointRtree> &Rtree_list, short transi_level){
    if (node_level == transi_level+1){
        Rtree_list.push_back(Rtree);
        return;
    }

    if (LUNode){
        LUNode->locRtree(Rtree_list, transi_level);
    }

    if (RUNode){
        RUNode->locRtree(Rtree_list, transi_level);
    }

    if (LBNode){
        LBNode->locRtree(Rtree_list, transi_level);
    }

    if (RBNode){
        RBNode->locRtree(Rtree_list, transi_level);
    }

    return;
}


// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// **********************************************    linestring node functions    *************************************************************** //
// ---------------------------------------------------------------------------------------------------------------------------------------------- //
linestringNode::linestringNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, linestringNode *_parent){
    this->node_type = _node_type;
    this->node_x = _node_x;
    this->node_y = _node_y;
    this->node_w = _node_w;
    this->node_level = _node_level;
    this->node_fcount = 0;
    this->node_scount = 0;
    this->node_flen = 0;
    this->parent = _parent;
    this->LUNode = nullptr;
    this->RUNode = nullptr;
    this->LBNode = nullptr;
    this->RBNode = nullptr;
}


linestringNode::~linestringNode(){
    delete LUNode;
    delete RUNode;
    delete LBNode;
    delete RBNode;
    delete parent;
    LUNode = nullptr;
    RUNode = nullptr;
    LBNode = nullptr;
    RBNode = nullptr;
    parent = nullptr;
}


// insert TMBB of a linestring feature
linestringNode *linestringNode::insertTMBB(double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax, int seg_count, double line_len, short transi_level){
    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= transi_level)){
        if (!LUNode){
            LUNode = new linestringNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        LUNode->node_fcount ++;
        LUNode->node_scount += seg_count;
        LUNode->node_flen += line_len;
        linestringNode *seek = LUNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, seg_count, line_len, transi_level);
        return seek;
    }

    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= transi_level)){
        if (!RUNode){
            RUNode = new linestringNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level+1, this);
        }
        RUNode->node_fcount ++;
        RUNode->node_scount += seg_count;
        RUNode->node_flen += line_len;
        linestringNode *seek = RUNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, seg_count, line_len, transi_level);
        return seek;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= transi_level)){
        if (!LBNode){
            LBNode = new linestringNode(LB, node_x, node_y, node_w/2, node_level+1, this);
        }
        LBNode->node_fcount ++;
        LBNode->node_scount += seg_count;
        LBNode->node_flen += line_len;
        linestringNode *seek = LBNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, seg_count, line_len, transi_level);
        return seek;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= transi_level)){
        if (!RBNode){
            RBNode = new linestringNode(RB, node_x + node_w/2, node_y, node_w/2, node_level+1, this);
        }
        RBNode->node_fcount ++;
        RBNode->node_scount += seg_count;
        RBNode->node_flen += line_len;
        linestringNode *seek = RBNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, seg_count, line_len, transi_level);
        return seek;
    }
    else {
        return this;
    }
}


// insert the segement of a linestring feature
void linestringNode::insertSeg(double seg_xMin, double seg_yMin, double seg_xMax, double seg_yMax, string seg_type, double line_len, short transi_level){    
    if (node_level == transi_level + 8){
        return;
    }

    // when the node contain the segement
    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, seg_xMin, seg_yMin, seg_xMax, seg_yMax)){
        if (!LUNode){
            LUNode = new linestringNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        LUNode->node_scount++;
        LUNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type, line_len, transi_level);
        return;
    }
    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, seg_xMin, seg_yMin, seg_xMax, seg_yMax)){
        if (!RUNode){
            RUNode = new linestringNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level + 1, this);
        }
        RUNode->node_scount++;
        RUNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type, line_len, transi_level);
        return;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, seg_xMin, seg_yMin, seg_xMax, seg_yMax)){
        if (!LBNode){
            LBNode = new linestringNode(LB, node_x, node_y, node_w/2, node_level + 1, this);
        }
        LBNode->node_scount++;
        LBNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type, line_len, transi_level);
        return;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, seg_xMin, seg_yMin, seg_xMax, seg_yMax)){
        if (!RBNode){
            RBNode = new linestringNode(RB, node_x + node_w/2, node_y, node_w/2, node_level + 1, this);
        }
        RBNode->node_scount++;
        RBNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type, line_len, transi_level);
        return;
    }
    // when the node overlap with the segement
    else{
        if (IfInterctSeg(node_x, node_y + node_w/2, node_w/2, seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type)){
            if (!LUNode){
                LUNode = new linestringNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
            }
            LUNode->node_scount++;
            LUNode->node_flen += line_len;
            LUNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type, line_len, transi_level);
        }
        if (IfInterctSeg(node_x + node_w/2, node_y+node_w/2, node_w/2, seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type)){
            if (!RUNode){
                RUNode = new linestringNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level + 1, this);
            }
            RUNode->node_scount++;
            RUNode->node_flen += line_len;
            RUNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type, line_len, transi_level);
        }
        if (IfInterctSeg(node_x, node_y, node_w/2, seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type)){
            if (!LBNode){
                LBNode = new linestringNode(LB, node_x, node_y, node_w/2, node_level + 1, this);
            }
            LBNode->node_scount++;
            LBNode->node_flen += line_len;
            LBNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type, line_len, transi_level);
        }
        if (IfInterctSeg(node_x + node_w/2, node_y, node_w/2, seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type)){
            if (!RBNode){
                RBNode = new linestringNode(RB, node_x + node_w/2, node_y, node_w/2, node_level + 1, this);
            }
            RBNode->node_scount++;
            RBNode->node_flen += line_len;
            RBNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, seg_type, line_len, transi_level);
        }
    }
    return;
}


void linestringNode::insertFeature(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short transi_level, BOX_PAIR line){
    if (node_level == transi_level + 1){
        Rtree.insert(line);
    }
    
    if (node_level == transi_level + 8){
        return;
    }


    // LU sub-region
    if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (LUNode){
            LUNode->node_fcount++;
            LUNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, transi_level, line);
        }
    }
    // RU sub-region
    if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (RUNode){
            RUNode->node_fcount++;
            RUNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, transi_level, line);
        }
    }
    // LB sub-region
    if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (LBNode){
            LBNode->node_fcount++;
            LBNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, transi_level, line);
        }
    }
    // RB sub-region
    if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (RBNode){
            RBNode->node_fcount++;
            RBNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax,transi_level, line);
        }
    }
    return; 
}


// loc the R-tree list of a node
void linestringNode::locRtree(vector<linestringRtree> &Rtree_list, short transi_level){
    if (node_level == transi_level+1){
        Rtree_list.push_back(Rtree);
        return;
    }

    if (LUNode){
        LUNode->locRtree(Rtree_list, transi_level);
    }

    if (RUNode){
        RUNode->locRtree(Rtree_list, transi_level);
    }

    if (LBNode){
        LBNode->locRtree(Rtree_list, transi_level);
    }

    if (RBNode){
        RBNode->locRtree(Rtree_list, transi_level);
    }

    return;
}





// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// ***************************************************    polygon node functions    *************************************************************** //
// ---------------------------------------------------------------------------------------------------------------------------------------------- //
polygonNode::polygonNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, polygonNode *_parent){
    this->node_type = _node_type;
    this->topo_type = 'i';
    this->node_x = _node_x;
    this->node_y = _node_y;
    this->node_w = _node_w;
    this->node_level = _node_level;
    this->node_fcount = 0;
    this->node_flen = 0;
    this->node_farea = 0;
    this->parent = _parent;
    this->LUNode = nullptr;
    this->RUNode = nullptr;
    this->LBNode = nullptr;
    this->RBNode = nullptr;
}


polygonNode::~polygonNode(){
    delete LUNode;
    delete RUNode;
    delete LBNode;
    delete RBNode;
    delete parent;
    LUNode = nullptr;
    RUNode = nullptr;
    LBNode = nullptr;
    RBNode = nullptr;
    parent = nullptr;
}


// insert TMBB of a polygon feature
polygonNode * polygonNode::insertTMBB(double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax, int edge_count, double obj_len, double obj_area, short transi_level){
    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= transi_level)){
        if (!LUNode){
            LUNode = new polygonNode(LU, node_x, node_y + node_w/2, node_w/2, node_level+1, this);
        }
        LUNode->node_fcount ++;
        LUNode->node_flen += obj_len;
        LUNode->node_farea += obj_area;
        polygonNode *seek = LUNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, edge_count, obj_len, obj_area, transi_level);
        return seek;
    }
    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= transi_level)){
        if (!RUNode){
            RUNode = new polygonNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level+1, this);
        }
        RUNode->node_fcount ++;
        RUNode->node_flen += obj_len;
        RUNode->node_farea += obj_area;
        polygonNode *seek = RUNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, edge_count, obj_len, obj_area, transi_level);
        return seek;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= transi_level)){
        if (!LBNode){
            LBNode = new polygonNode(LB, node_x, node_y, node_w/2, node_level+1, this);
        }
        LBNode->node_fcount ++;
        LBNode->node_flen += obj_len;
        LBNode->node_farea += obj_area;
        polygonNode *seek = LBNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, edge_count, obj_len, obj_area, transi_level);
        return seek;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= transi_level)){
        if (!RBNode){
            RBNode = new polygonNode(RB, node_x+node_w/2, node_y, node_w/2, node_level+1, this);
        }
        RBNode->node_fcount ++;
        RBNode->node_flen += obj_len;
        RBNode->node_farea += obj_area;
        polygonNode *seek = RBNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, edge_count, obj_len, obj_area, transi_level);
        return seek;
    }
    else {
        return this;
    }
}


// insert the edge of a polygon feature
void polygonNode::insertEdge(double edge_xMin, double edge_yMin, double edge_xMax, double edge_yMax, string edge_type, double obj_len, short transi_level){
    if (node_level == transi_level + 8){
        return;
    }

    // when the node contain the edge
    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax)){
        if (!LUNode){
            LUNode = new polygonNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        
        LUNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, obj_len, transi_level);
        return;
    }
    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax)){
        if (!RUNode){
            RUNode = new polygonNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level+1, this);
        }
        
        RUNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, obj_len, transi_level);
        return;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax)){
        if (!LBNode){
            LBNode = new polygonNode(LB, node_x, node_y, node_w/2, node_level+1, this);
        }
        
        LBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, obj_len, transi_level);
        return;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax)){
        if (!RBNode){
            RBNode = new polygonNode(RB, node_x+node_w/2, node_y, node_w/2, node_level+1, this);
        }
        
        RBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, obj_len, transi_level);
        return;
    }
    // when the node overlap with the edge
    else{
        if (IfInterctSeg(node_x, node_y + node_w/2, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type)){
            if (!LUNode){
                LUNode = new polygonNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
            }
            
            LUNode->node_flen += obj_len;
            LUNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, obj_len, transi_level);
        }
        if (IfInterctSeg(node_x + node_w/2, node_y+node_w/2, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type)){
            if (!RUNode){
                RUNode = new polygonNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level + 1, this);
            }
            
            RUNode->node_flen += obj_len;
            RUNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, obj_len, transi_level);
        }
        if (IfInterctSeg(node_x, node_y, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type)){
            if (!LBNode){
                LBNode = new polygonNode(LB, node_x, node_y, node_w/2, node_level + 1, this);
            }

            LBNode->node_flen += obj_len;
            LBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, obj_len, transi_level);
        }
        if (IfInterctSeg(node_x + node_w/2, node_y, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type)){
            if (!RBNode){
                RBNode = new polygonNode(RB, node_x + node_w/2, node_y, node_w/2, node_level + 1, this);
            }
            
            RBNode->node_flen += obj_len;
            RBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, obj_len, transi_level);
        }
    }
    return;
}


// insert the polygon feature
void polygonNode::insertFeature(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short transi_level, double ring_x[], double ring_y[], int len, double obj_area, BOX_PAIR poly){
    if (node_level == transi_level + 1){
        Rtree.insert(poly);
    }
    
    if (node_level == transi_level + 8 || topo_type == 'w'){
        return;
    }

    // LU sub-region
    if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LUNode){
            if (IfwithinBox(node_x+node_w/4, node_y+3*node_w/4, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (IfwithinPolygon(node_x+node_w/4, node_y+3*node_w/4, ring_x, ring_y, len)){
                    LUNode = new polygonNode(LU, node_x, node_y+node_w/2, node_w/2, node_level+1, this);
                    LUNode->topo_type = 'w';
                    LUNode->node_fcount++;
                }
            }
        }
        else{
            LUNode->node_fcount++;
            LUNode->node_farea += obj_area;
            LUNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, transi_level, ring_x, ring_y, len, obj_area, poly);
        }
    }
    // RU sub-region
    if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RUNode){
            if (IfwithinBox(node_x+3*node_w/4, node_y+3*node_w/4, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (IfwithinPolygon(node_x+3*node_w/4, node_y+3*node_w/4, ring_x, ring_y, len)){
                    RUNode = new polygonNode(RU, node_x+node_w/2, node_y+node_w/2, node_w/2, node_level+1, this);
                    RUNode->topo_type = 'w';
                    RUNode->node_fcount++;
                }
            }
        }
        else{
            RUNode->node_fcount++;
            RUNode->node_farea += obj_area;
            RUNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, transi_level, ring_x, ring_y, len, obj_area, poly);
        } 
    }
    // LB sub-region
    if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LBNode){
            if (IfwithinBox(node_x+node_w/4, node_y+node_w/4, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (IfwithinPolygon(node_x+node_w/4, node_y+node_w/4, ring_x, ring_y, len)){
                    LBNode = new polygonNode(LB, node_x, node_y, node_w/2, node_level+1, this);
                    LBNode->topo_type = 'w';
                    LBNode->node_fcount++;
                }  
            }
        }
        else{
            LBNode->node_fcount++;
            LBNode->node_farea += obj_area;
            LBNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, transi_level, ring_x, ring_y, len, obj_area, poly);
        }
    }
    // RB sub-region
    if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RBNode){
            if (IfwithinBox(node_x+3*node_w/4, node_y+node_w/4, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (IfwithinPolygon(node_x+3*node_w/4, node_y+node_w/4, ring_x, ring_y, len)){
                    RBNode = new polygonNode(RB, node_x+node_w/2, node_y, node_w/2, node_level+1, this);
                    RBNode->topo_type = 'w';
                    RBNode->node_fcount++;
                } 
            }
        }
        else{
            RBNode->node_fcount++;
            RBNode->node_farea += obj_area;
            RBNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax,transi_level, ring_x, ring_y, len, obj_area, poly);
        }
    }
    return; 
}


// loc the R-tree list of a node
void polygonNode::locRtree(vector<polygonRtree> &Rtree_list, short transi_level){
    if (node_level == transi_level+1 && topo_type == 'i'){
        Rtree_list.push_back(Rtree);
        return;
    }

    if (LUNode){
        LUNode->locRtree(Rtree_list, transi_level);
    }

    if (RUNode){
        RUNode->locRtree(Rtree_list, transi_level);
    }

    if (LBNode){
        LBNode->locRtree(Rtree_list, transi_level);
    }

    if (RBNode){
        RBNode->locRtree(Rtree_list, transi_level);
    }

    return;
}


void polygonNode::searchTile(int z, int x, int y, int count, int transi_level, vector<int> &tile_id){
    if (node_level == transi_level + 8 || node_fcount < count){
        if (z != 0){
            return;
        }
    }

    if (node_fcount >= count and node_fcount <= count + 100){
        tile_id.push_back(z);
        tile_id.push_back(x);
        tile_id.push_back(y);
        return;
    }

    if (LUNode){
        LUNode->searchTile(z+1, 2*x, 2*y, count, transi_level, tile_id);
    }

    if (RUNode){
        RUNode->searchTile(z+1, 2*x+1, 2*y, count, transi_level, tile_id);
    }

    if (LBNode){
        LBNode->searchTile(z+1, 2*x, 2*y+1, count, transi_level, tile_id);
    }

    if (RBNode){
        RBNode->searchTile(z+1, 2*x+1, 2*y+1, count, transi_level, tile_id);
    }

    return;
}

