#include <chrono>
#include <filesystem>
#include <sys/stat.h>
#include <random>
#include <numeric>
#include <mpi.h>
#include "buildIndex.hpp"
#include "pointRender.hpp"
#include "linestringRender.hpp"
#include "polygonRender.hpp"

namespace fs = filesystem;
#define L 20037508.3427892

typedef bgm::polygon<POINT> polygon_type;
typedef bgm::linestring<POINT> line_type;
typedef bgi::rtree<pair<Box, polygon_type>, bgi::quadratic<8>> polyRtree;
typedef bgi::rtree<pair<Box, line_type>, bgi::quadratic<8>> lineRtree;

struct line_within_box{
    const Box& box;

    line_within_box(const Box& box) : box(box) {}

    bool operator()(const BOX_PAIR& value) const {
        return bg::within(value.first, box);
    }
};


struct poly_within_box{
    const Box& box;

    poly_within_box(const Box& box) : box(box) {}

    bool operator()(const BOX_PAIR& poly) const {
        return bg::within(poly.first, box);
    }
};


struct KeyValue {
    int id;
    int count;
    double total_time;
};


void mergeTime(int process_num, map<int, vector<double>> time_map, map<int, double> &speed_map){
    vector<KeyValue> part_vec;
    for (const auto& pair: time_map) {
        vector<double> time_vec = pair.second;
        if (time_vec.size() != 0){
            double total_time = accumulate(time_vec.begin(), time_vec.end(), 0.0);
            part_vec.push_back({pair.first, time_vec.size(), total_time});
        }
        else{
            part_vec.push_back({pair.first, 0, 0});
        }
    }
    int local_size = part_vec.size();
    int global_size = 0;

    MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    vector<KeyValue> global_vec(global_size);

    MPI_Datatype mpi_type;
    int blocklengths[3] = {1, 1, 1};
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};

    displacements[0] = offsetof(KeyValue, id);
    displacements[1] = offsetof(KeyValue, count);
    displacements[2] = offsetof(KeyValue, total_time);

    MPI_Type_create_struct(3, blocklengths, displacements, types, &mpi_type);
    MPI_Type_commit(&mpi_type);

    MPI_Allgather(part_vec.data(), local_size, mpi_type, global_vec.data(), local_size, mpi_type, MPI_COMM_WORLD);


    map<int, int> total_count_map;
    map<int, double> total_time_map;
    for (auto kv: global_vec) {
        total_count_map[kv.id] += kv.count;
        total_time_map[kv.id] += kv.total_time;
    }

    
    for (int i = 0; i < 6; i++){
        int count = total_count_map[i];
        if (count != 0){
            speed_map[i] = (count * process_num) / total_time_map[i];
        }
        else{
            speed_map[i] = 0;
        }
    }

    MPI_Type_free(&mpi_type);
}


line_type createLine(vector<double> line_coords){
    line_type line;
    
    for (int i = 0; i < line_coords.size(); i+=2){
        line.push_back(POINT(line_coords[i], line_coords[i+1]));
    }

    return line;
}


polygon_type createPolygon(vector<double> poly_coords){
    polygon_type polygon;
    
    for (int i = 0; i < poly_coords.size(); i+=2){
        bg::append(polygon.outer(), POINT(poly_coords[i], poly_coords[i+1]));
    }

    return polygon;
}




lineRtree realRtree_line(vector<BOX_PAIR> lines){
    lineRtree line_rtree;

    for (int i = 0; i < lines.size(); i++){
        line_type line = createLine(lines[i].second);
        Box line_box;
        bg::envelope(line, line_box);
        line_rtree.insert(make_pair(line_box, line));
    }
    return line_rtree;
}


polyRtree realRtree_poly(vector<BOX_PAIR> polys){
    polyRtree poly_rtree;

    for (int i = 0; i < polys.size(); i++){
        polygon_type polygon = createPolygon(polys[i].second);
        Box poly_box;
        bg::envelope(polygon, poly_box);
        poly_rtree.insert(make_pair(poly_box, polygon));
    }

    // polygon_type f_polygon = createPolygon(polys[0].second);
    // Box f_poly_box;
    // bg::envelope(f_polygon, f_poly_box);
    // poly_rtree.insert(make_pair(f_poly_box, f_polygon));
    
    // for (int i = 1; i < polys.size(); i++){
    //     polygon_type s_polygon = createPolygon(polys[i].second);
    //     if (!bg::equals(f_polygon, s_polygon)){
    //         Box poly_box;
    //         bg::envelope(s_polygon, poly_box);
    //         poly_rtree.insert(make_pair(poly_box, s_polygon));

    //         f_polygon = s_polygon;
    //     }
    // }

    return poly_rtree;
}


bool getShpFiles(string data_path, vector<string> &shp_files){
    // check the shapefile folder if exists
    struct stat buffer;
    if (stat(data_path.c_str(), &buffer) != 0){
        cout << "【ERROR】 please check the shapefile folder whether exists !" << endl;
        cout << "【ERROR】 " << data_path << endl;
        return false;
    }
    
    // find the .shp file from the folder
    for (const auto& entry: fs::recursive_directory_iterator(data_path)){
        if (entry.is_regular_file() && entry.path().extension() == ".shp") {
            shp_files.push_back(entry.path().string());
        }
    }

    if (shp_files.size()){
        return true;
    }
    else{
        cout << "【ERROR】 please check the shapefile folder whether contains .shp file !" << endl;
        cout << "【ERROR】 " << data_path << endl;
        return false;
    }
}


map<int, vector<int>> calTileBox(vector<string> shp_files, short min_level, short max_level){
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
    
    
    coordTran(data_xmin, data_ymin);
    coordTran(data_xmax, data_ymax);

    map<int, vector<int>> tile_box_map;
    for (int level = min_level; level <= max_level; level++){
        double tile_Rz = L / pow(2, level-1);
        int tile_xMin = floor((data_xmin + L) / tile_Rz);
        int tile_xMax = floor((data_xmax + L) / tile_Rz);
        int tile_yMin = floor((L - data_ymax) / tile_Rz);
        int tile_yMax = floor((L - data_ymin) / tile_Rz);

        tile_box_map[level] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    }

    return tile_box_map;
}


double calANN(int z, pointNode *tile_node, short transi_level){
    int feature_count = tile_node->node_fcount;

    if (feature_count < 2 or feature_count > 2500){
        return 0;
    }

    double NNI = 0;
    
    double tile_Rz = L / pow(2, z-1);

    vector<pointRtree> Rtree_list;
    tile_node->locRtree(Rtree_list, transi_level);

    pointRtree tile_rtree;
    for (auto rtree: Rtree_list){
        for (const auto& item : rtree) {
            tile_rtree.insert(item);
        }
    }

    for (const auto& obj: tile_rtree){
        vector<POINT> near_pt;
        tile_rtree.query(bgi::nearest(obj, 2), back_inserter(near_pt));
        NNI += bg::distance(obj, near_pt[0]) / (1000*feature_count);
    }

    double ANN = 2 * NNI * sqrt(feature_count / (tile_Rz*tile_Rz/1000));

    return ANN;
}


double calANN(int z, int x, int y, linestringNode *tile_node, short transi_level){
    int feature_count = tile_node->node_fcount;

    if (feature_count < 2 or feature_count > 5000){
        return 0;
    }

    double NNI = 0;
    
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;

    vector<linestringRtree> Rtree_list;
    tile_node->locRtree(Rtree_list, transi_level);

    pointRtree tile_rtree;
    for (auto rtree: Rtree_list){
        for (const auto& line : rtree) {
            POINT cent_point;
            bg::centroid(createLine(line.second), cent_point);
            tile_rtree.insert(cent_point);
        }
    }

    for (const auto& obj: tile_rtree){
        vector<POINT> near_pt;
        tile_rtree.query(bgi::nearest(obj, 2), back_inserter(near_pt));
        NNI += bg::distance(obj, near_pt[0]) / (1000*feature_count);
    }

    double ANN = 2 * NNI * sqrt(feature_count / (tile_Rz*tile_Rz/1000));

    return ANN;
}


double calANN(int z, int x, int y, polygonNode *tile_node, short transi_level){
    int feature_count = tile_node->node_fcount;

    if (feature_count < 2 or feature_count > 2000){
        return 0;
    }

    double NNI = 0;
    
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;

    vector<polygonRtree> Rtree_list;
    tile_node->locRtree(Rtree_list, transi_level);

    pointRtree tile_rtree;
    for (auto rtree: Rtree_list){
        for (const auto& poly : rtree) {
            POINT cent_point;
            bg::centroid(createPolygon(poly.second), cent_point);
            tile_rtree.insert(cent_point);
        }
    }

    for (const auto& obj: tile_rtree){
        vector<POINT> near_pt;
        tile_rtree.query(bgi::nearest(obj, 2), back_inserter(near_pt));
        NNI += bg::distance(obj, near_pt[0]) / (1000*feature_count);
    }

    double ANN = 2 * NNI * sqrt(feature_count / (tile_Rz*tile_Rz/1000));

    return ANN;
}


void calSpeed(int cpu_rank, int process_num, vector<string> shp_files, pointNode *RNode, double min_ann, double max_ann, short min_level, short max_level){
    map<int, vector<double>> DiPG_time_map, DaPG_time_map;
    for (int i = 0; i < 6; i++){
        DiPG_time_map[i] = vector<double>();
        DaPG_time_map[i] = vector<double>();
    }
    int total_num = 2000;
    int tile_num = 0;

    map<int, vector<int>> tile_box_map = calTileBox(shp_files, min_level, max_level);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> z_dis(min_level, max_level);
    
    while (1){
        int z = z_dis(gen);
        vector<int> tile_box = tile_box_map[z];
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);
        int x = x_dis(gen);
        int y = y_dis(gen);

        pointNode *tile_node = locTileNode(z, x, y, RNode);
        if (tile_node){
            double f_ann = calANN(z, tile_node, max_level);
            if (abs(f_ann - 0) <= 1e-9){
                continue;
            }

            double DiPG_time = DiPGPlot(z, x, y, 2, RNode);   
            double DaPG_time = DaPGPlot(z, x, y, 2, RNode, max_level);

            if (f_ann <= min_ann){
                DiPG_time_map[0].push_back(DiPG_time);
                DaPG_time_map[0].push_back(DaPG_time);
            }
            else if (f_ann > max_ann){
                DiPG_time_map[5].push_back(DiPG_time);
                DaPG_time_map[5].push_back(DaPG_time);
            }
            else{
                int ix = static_cast<int>(ceil(1000*(f_ann-min_ann) / 5.0));
                DiPG_time_map[ix].push_back(DiPG_time);
                DaPG_time_map[ix].push_back(DaPG_time);
            }

            tile_num++;
        }

        if (tile_num >= total_num / process_num){
            if (DiPG_time_map[0].size() and DiPG_time_map[1].size() and DiPG_time_map[2].size()){
                if (DiPG_time_map[3].size() and DiPG_time_map[4].size() and DiPG_time_map[5].size()){
                    break;
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    map<int, double> DiPG_speed_map, DaPG_speed_map;
    mergeTime(process_num, DiPG_time_map, DiPG_speed_map);
    mergeTime(process_num, DaPG_time_map, DaPG_speed_map);


    if (cpu_rank == 0){
        for (int i = 0; i < 6; i++){
            cout << DiPG_speed_map[i] << " " << DaPG_speed_map[i] << endl;
        }
    }
}


void calSpeed(int cpu_rank, int process_num, vector<string> shp_files, linestringNode *RNode, double min_ann, double max_ann, short min_level, short max_level){
    map<int, vector<double>> DiPG_time_map, DaPG_time_map;
    for (int i = 0; i < 6; i++){
        DiPG_time_map[i] = vector<double>();
        DaPG_time_map[i] = vector<double>();
    }
    int total_num = 2000;
    int tile_num = 0;

    map<int, vector<int>> tile_box_map = calTileBox(shp_files, min_level, max_level);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> z_dis(min_level, max_level);
    
    while (1){
        int z = z_dis(gen);
        vector<int> tile_box = tile_box_map[z];
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);
        int x = x_dis(gen);
        int y = y_dis(gen);

        linestringNode *tile_node = locTileNode(z, x, y, RNode);
        if (tile_node){
            double f_ann = calANN(z, x, y, tile_node, max_level);
            if (abs(f_ann - 0) <= 1e-9){
                continue;
            }

            double DiPG_time = DiPGPlot(z, x, y, 1, RNode);   
            double DaPG_time = DaPGPlot(z, x, y, 2, RNode, max_level);

            if (f_ann <= min_ann){
                DiPG_time_map[0].push_back(DiPG_time);
                DaPG_time_map[0].push_back(DaPG_time);
            }
            else if (f_ann > max_ann){
                DiPG_time_map[5].push_back(DiPG_time);
                DaPG_time_map[5].push_back(DaPG_time);
            }
            else{
                int ix = static_cast<int>(ceil(1000*(f_ann-min_ann) / 5.0));
                DiPG_time_map[ix].push_back(DiPG_time);
                DaPG_time_map[ix].push_back(DaPG_time);
            }

            tile_num++;
        }

        if (tile_num >= total_num / process_num){
            if (DiPG_time_map[0].size() and DiPG_time_map[1].size() and DiPG_time_map[2].size()){
                if (DiPG_time_map[3].size() and DiPG_time_map[4].size() and DiPG_time_map[5].size()){
                    break;
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    map<int, double> DiPG_speed_map, DaPG_speed_map;
    mergeTime(process_num, DiPG_time_map, DiPG_speed_map);
    mergeTime(process_num, DaPG_time_map, DaPG_speed_map);


    if (cpu_rank == 0){
        for (int i = 0; i < 6; i++){
            cout << DiPG_speed_map[i] << " " << DaPG_speed_map[i] << endl;
        }
    }
}


void calSpeed(int cpu_rank, int process_num, vector<string> shp_files, polygonNode *RNode, double min_ann, double max_ann, short min_level, short max_level){
    map<int, vector<double>> DiPG_time_map, DaPG_time_map;
    for (int i = 0; i < 6; i++){
        DiPG_time_map[i] = vector<double>();
        DaPG_time_map[i] = vector<double>();
    }
    int total_num = 2000;
    int tile_num = 0;

    map<int, vector<int>> tile_box_map = calTileBox(shp_files, min_level, max_level);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> z_dis(min_level, max_level);
    
    while (1){
        int z = z_dis(gen);
        vector<int> tile_box = tile_box_map[z];
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);
        int x = x_dis(gen);
        int y = y_dis(gen);

        polygonNode *tile_node = locTileNode(z, x, y, RNode);
        
        if (tile_node && tile_node->topo_type == 'i'){
            double f_ann = calANN(z, x, y, tile_node, max_level);
            if (abs(f_ann - 0) <= 1e-9){
                continue;
            }
            
            double DiPG_time = DiPGPlot(z, x, y, 1, RNode);   
            double DaPG_time = DaPGPlot(z, x, y, 2, RNode, max_level);

            if (f_ann <= min_ann){
                DiPG_time_map[0].push_back(DiPG_time);
                DaPG_time_map[0].push_back(DaPG_time);
            }
            else if (f_ann > max_ann){
                DiPG_time_map[5].push_back(DiPG_time);
                DaPG_time_map[5].push_back(DaPG_time);
            }
            else{
                int ix = static_cast<int>(ceil(1000*(f_ann-min_ann) / 5.0));
                DiPG_time_map[ix].push_back(DiPG_time);
                DaPG_time_map[ix].push_back(DaPG_time);
            }

            tile_num++;
        }

        if (tile_num >= total_num / process_num){
            if (DiPG_time_map[0].size() and DiPG_time_map[1].size() and DiPG_time_map[2].size()){
                if (DiPG_time_map[3].size() and DiPG_time_map[4].size() and DiPG_time_map[5].size()){
                    break;
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    map<int, double> DiPG_speed_map, DaPG_speed_map;
    mergeTime(process_num, DiPG_time_map, DiPG_speed_map);
    mergeTime(process_num, DaPG_time_map, DaPG_speed_map);


    if (cpu_rank == 0){
        for (int i = 0; i < 6; i++){
            cout << DiPG_speed_map[i] << " " << DaPG_speed_map[i] << endl;
        }
    }
}





int main(int argc, char **argv){
    int cpu_rank, process_num; 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &cpu_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);
    
    string data_path = argv[1];
    string shp_type = argv[2];
    double min_ann = stod(argv[3]);
    double max_ann = stod(argv[4]);
    short min_level = stoi(argv[5]);
    short max_level = stoi(argv[6]);

    vector<string> shp_files;

    if (getShpFiles(data_path, shp_files)){
        if (shp_type == "point"){
            pointNode *pointRNode = pointIndex(shp_files, "4326", max_level);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of DiPG and DaPG on diffirent feature ANN interval
            calSpeed(cpu_rank, process_num, shp_files, pointRNode, min_ann, max_ann, min_level, max_level);
        }
        else if (shp_type == "linestring"){
            linestringNode *linestringRNode = linestringIndex(shp_files, "4326", max_level);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of DiPG and DaPG on diffirent feature ANN interval
            calSpeed(cpu_rank, process_num, shp_files, linestringRNode, min_ann, max_ann, min_level, max_level);
        }
        else if (shp_type == "polygon"){
            polygonNode *polygonRNode = polygonIndex(shp_files, "4326", max_level);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of DiPG and DaPG on diffirent feature ANN interval
            calSpeed(cpu_rank, process_num, shp_files, polygonRNode, min_ann, max_ann, min_level, max_level);
        }
    }

    MPI_Finalize();
    return 0;
}