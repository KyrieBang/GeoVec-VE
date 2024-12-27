#include "buildIndex.hpp"
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <iterator>
#include <limits>

#define pi 3.14159265358979323
#define L 20037508.3427892


double calArea(double *x_list, double *y_list, int coord_num){
    double area = 0.0;
    for (int i = 0; i < coord_num-1; ++i) {
        area += x_list[i] * y_list[i+1] - y_list[i] * x_list[i+1];
    }
    area += x_list[coord_num-1] * y_list[0] - y_list[coord_num-1] * x_list[0];

    return abs(area) / (2.0 * 1000000);
}

// coordinate transformation (4326->3857)
void coordTran(double &x ,double &y){
    x = x * 20037508.34 / 180;
    y = (log(tan(((90 + y) * pi) / 360)) / (pi / 180)) * 20037508.34 / 180;
}


// calculate the spatial scale of the dataset
short calTransiLevel(vector<string> shp_files, string data_espg){
    Shapefile shp;
    short data_level = 0;
    double data_xmin = numeric_limits<double>::infinity();
    double data_xmax = -numeric_limits<double>::infinity();
    double data_ymin = numeric_limits<double>::infinity();
    double data_ymax = -numeric_limits<double>::infinity();

    for (auto shp_path: shp_files){
        GeoData data = shp.Read(shp_path);

        double shp_xmin = data.extent.xmin;
        double shp_xmax = data.extent.xmax;
        double shp_ymin = data.extent.ymin;
        double shp_ymax = data.extent.ymax;

        data_xmin = min(data_xmin, shp_xmin);
        data_xmax = max(data_xmax, shp_xmax);
        data_ymin = min(data_ymin, shp_ymin);
        data_ymax = max(data_ymax, shp_ymax);
    }
    
    
    if (data_espg == "4326"){
        coordTran(data_xmin, data_ymin);
        coordTran(data_xmax, data_ymax);
    }


    double dlt_x = data_xmax - data_xmin;
    double dlt_y = data_ymax - data_ymin;
    double dlt_max = max(dlt_x, dlt_y);
    double span = 2 * L;
    while (dlt_max < span){
        span = span / 2;
        data_level++;
    } 
    return data_level;
}


// build index of point dataset
pointNode* pointIndex(vector<string> shp_files, string data_espg, short transi_level){
    // init
    Shapefile shp;
    int data_count = 0;
    
    struct timeval	t1, t2;
	gettimeofday(&t1, NULL);

    // create root node ptr 
    pointNode *pointRNode = new pointNode(ROOT, -L, -L, 2 * L, 0, nullptr);
    
    // read shapfile
    for (auto shp_path: shp_files){
        GeoData data = shp.Read(shp_path);
        
        // traverse point feature
        PointIterator ptIter(data);
        do{
            double x = ptIter.X();
            double y = ptIter.Y();
            
            if (data_espg == "4326"){
                coordTran(x, y);
            }
            // insert point from root node
            pointRNode->insertObj(x, y, transi_level);
            data_count++;
        } while (ptIter.NextFeature() >= 0);
    }

    gettimeofday(&t2, NULL);
    float time_use = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;
    
    return pointRNode;
}



// build index of linestring dataset
linestringNode* linestringIndex(vector<string> shp_files, string data_espg, short transi_level){
    // init
    Shapefile shp;
    int feature_count = 0;
    int segment_count = 0;
    
    struct timeval	t1, t2;
	gettimeofday(&t1, NULL);

    // create root node ptr
    linestringNode *lineRNode = new linestringNode(ROOT, -L, -L, 2 * L, 0, nullptr);

    // read shapefile
    for (auto shp_path: shp_files){
        GeoData data = shp.Read(shp_path);
        
        // traverse linestring feature
        LineStringIterator lineIter(data);
        do{
            int nodes_num = lineIter.NodeCount();
            double line_x[nodes_num];
            double line_y[nodes_num];
            XY_VEC line_cord;
            double *nodes = lineIter.Nodes();
            
            for (int i = 0; i < nodes_num; i++){
                double xi = lineIter.X(nodes, i);
                double yi = lineIter.Y(nodes, i);
                if (data_espg == "4326"){
                    coordTran(xi, yi);
                }
                line_x[i] = xi;
                line_y[i] = yi;
                line_cord.push_back(xi);
                line_cord.push_back(yi);
            }
            double x_min = *min_element(line_x, line_x+nodes_num);
            double x_max = *max_element(line_x, line_x+nodes_num);
            double y_min = *min_element(line_y, line_y+nodes_num);
            double y_max = *max_element(line_y, line_y+nodes_num);

            // cal the linestring length
            double line_len = 0;
            for (int k = 0; k < nodes_num-1; k++){
                line_len += (sqrt(pow((line_x[k+1]-line_x[k]), 2) + pow((line_y[k+1]-line_y[k]), 2)) / 1000.0);
            }

            // insert TMBB from root node
            linestringNode *MBNode = lineRNode->insertTMBB(x_min, y_min, x_max, y_max, (nodes_num-1), line_len, transi_level);
            // insert each segement from the minimum bounding node
            for (int j = 0; j < nodes_num - 1; j++){
                double x0 = line_x[j];
                double y0 = line_y[j];
                double x1 = line_x[j + 1];
                double y1 = line_y[j + 1];
                double seg_xMin = min(x0, x1);
                double seg_xMax = max(x0, x1);
                double seg_yMin = min(y0, y1);
                double seg_yMax = max(y0, y1);
                if (((x1 >= x0) && (y1 >= y0)) or ((x1 <= x0) && (y1 <= y0))){
                    MBNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, "0", line_len, transi_level);
                }
                else if (((x1 < x0) && (y1 > y0)) or ((x1 > x0) && (y1 < y0))){
                    MBNode->insertSeg(seg_xMin, seg_yMin, seg_xMax, seg_yMax, "1", line_len, transi_level);
                }
            }
            // insert the linestring feature from the minimum bounding node
            BOX_PAIR s_line = make_pair(Box(POINT(x_min, y_min), POINT(x_max, y_max)), line_cord); 
            MBNode->insertFeature(x_min, y_min, x_max, y_max, transi_level, s_line);
            
            segment_count += (nodes_num - 1);
            feature_count++;
        } while(lineIter.NextFeature() >= 0);
    }
    
    
    gettimeofday(&t2, NULL);
    float time_use = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;
    
    return lineRNode;
}




// build index of polygon dataset
polygonNode* polygonIndex(vector<string> shp_files, string data_espg, short transi_level){
    // init
    Shapefile shp;
    int feature_count = 0;
    int segment_count = 0;
    
    struct timeval	t1, t2;
	gettimeofday(&t1, NULL);

    // create root node
    polygonNode *polygonRNode = new polygonNode(ROOT, -L, -L, 2 * L, 0, nullptr);

    // read shapefile
    for (auto shp_path: shp_files){
        GeoData data = shp.Read(shp_path);
    
        // traverse polygon feature
        PolygonIterator polygonIter(data);
        do{
            // get the polygon feature
            for (int n = 0; n < polygonIter.PartCount(); n++){
                auto poly = polygonIter.Part(n);
                // get the ring of each polygon
                for (int k = 0; k < polygonIter.RingCount(poly); k++){
                    auto ring = polygonIter.Ring(poly, k);
                    auto nodes = polygonIter.Nodes(ring);
                    int nodes_num = polygonIter.NodeCount(ring);
                    double ring_x[nodes_num];
                    double ring_y[nodes_num];
                    XY_VEC ring_cord;
                    // for each node of each ring
                    for (int i = 0; i < nodes_num; i++){
                        double xi = polygonIter.X(nodes, i);
                        double yi = polygonIter.Y(nodes, i);
                        if (data_espg == "4326"){
                            coordTran(xi, yi);
                        }
                        ring_x[i] = xi;
                        ring_y[i] = yi;
                        ring_cord.push_back(xi);
                        ring_cord.push_back(yi);
                    }
                    double x_min = *min_element(ring_x, ring_x+nodes_num);
                    double x_max = *max_element(ring_x, ring_x+nodes_num);
                    double y_min = *min_element(ring_y, ring_y+nodes_num);
                    double y_max = *max_element(ring_y, ring_y+nodes_num);

                    // cal the ring length
                    double ring_len = 0;
                    for (int m = 0; m < nodes_num-1; m++){
                        ring_len += (sqrt(pow((ring_x[m+1]-ring_x[m]), 2) + pow((ring_y[m+1]-ring_y[m]), 2)) / 1000.0);
                    }
                    double ring_area = calArea(ring_x, ring_y, nodes_num-1);

                    if (ring_area <= 0){
                        continue;
                    }

                    // insert TMBB from root node
                    polygonNode *MBNode = polygonRNode->insertTMBB(x_min, y_min, x_max, y_max, (nodes_num-1), ring_len, ring_area, transi_level);
                    // insert each edge from the minimum bounding node
                    for (int j = 0; j < nodes_num - 1; j++){
                        double x0 = ring_x[j];
                        double y0 = ring_y[j];
                        double x1 = ring_x[j + 1];
                        double y1 = ring_y[j + 1];
                        double edge_xMin = min(x0, x1);
                        double edge_xMax = max(x0, x1);
                        double edge_yMin = min(y0, y1);
                        double edge_yMax = max(y0, y1);
                
                        if (((x1 >= x0) && (y1 >= y0)) or ((x1 <= x0) && (y1 <= y0))){
                            MBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, "0", ring_len, transi_level);
                        }
                        else if (((x1 < x0) && (y1 > y0)) or ((x1 > x0) && (y1 < y0))){
                            MBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, "1", ring_len, transi_level);
                        }
                    }
                    // insert the polygon feature from the minimum bounding node
                    BOX_PAIR s_poly = make_pair(Box(POINT(x_min, y_min), POINT(x_max, y_max)), ring_cord); 
                    MBNode->insertFeature(x_min, y_min, x_max, y_max, transi_level, ring_x, ring_y, nodes_num, ring_area, s_poly);
                    segment_count += (nodes_num-1);
                }
            }
        }while(polygonIter.NextFeature() >= 0);
    }
    
    gettimeofday(&t2, NULL);
    float time_use = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;
    
    return polygonRNode;
}


