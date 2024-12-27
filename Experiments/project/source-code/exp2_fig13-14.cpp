#include <filesystem>
#include <sys/stat.h>
#include <random>
#include <mpi.h>
#include "buildIndex.hpp"
#include "pointRender.hpp"
#include "linestringRender.hpp"
#include "polygonRender.hpp"

namespace fs = filesystem;
#define L 20037508.3427892

bool getShpFiles(string data_path, vector<string> &shp_files){
    // check the shapefile folder if exists
    struct stat buffer;
    if (stat(data_path.c_str(), &buffer) != 0){
        cout << "【ERROR】 please check the shapefile folder whether exists !" << endl;
        cout << "【ERROR】 " << data_path << endl;
        return false;
    }
    
    // find the .shp file from the folder
    for (const auto& entry: fs::directory_iterator(data_path)){
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


void calTimeSpeed(int cpu_rank, int process_num, vector<string> shp_files, pointNode *RNode, short level_thd, int count_thd){
    int total_num = 5000;
    int tile_num = 0;
    double total_time = 0;
    double max_time = 0;
    vector<double> part_time_list;
    vector<double> global_time_list(total_num);

    map<int, vector<int>> tile_box_map = calTileBox(shp_files, 0, 20);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> z_dis(0, 20);

    while (1){
        int z = z_dis(gen);
        vector<int> tile_box = tile_box_map[z];
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);
        int x = x_dis(gen);
        int y = y_dis(gen);

        
        pointNode *tile_node = locTileNode(z, x, y, RNode);
        if (tile_node){
            double tile_time = 0;

            if (tile_node->node_level <= level_thd){
                if (tile_node->node_fcount == 0){
                    continue;
                }

                if (tile_node->node_fcount < count_thd){
                    tile_time = DaPGPlot(z, x, y, 2, RNode, level_thd);
                }
                else{
                    tile_time = DiPGPlot(z, x, y, 2, RNode);
                }
            }
            else{
                auto s_time = chrono::high_resolution_clock::now();

                dataDrivenPlot(z, x, y, 255, 0, 0, 1, 2, RNode, level_thd);

                auto e_time = chrono::high_resolution_clock::now();

                chrono::duration<double> time_use = e_time - s_time;

                tile_time = time_use.count();
            }

            part_time_list.push_back(tile_time);

            total_time += tile_time;
            tile_num++;
        }

        if (tile_num == total_num / process_num){
            break;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Allgather(part_time_list.data(), total_num / process_num, MPI_DOUBLE, 
                  global_time_list.data(), total_num / process_num, MPI_DOUBLE, 
                  MPI_COMM_WORLD);
    
    
    if (cpu_rank == 0){
        for (auto tile_time: global_time_list){
            cout << tile_time << endl;
        }

        cout << max_time << " " << total_num / max_time << endl;
    }
}


void calTimeSpeed(int cpu_rank, int process_num, vector<string> shp_files, linestringNode *RNode, short level_thd, int count_thd, double len_thd){
    int total_num = 5000;
    int tile_num = 0;
    double total_time = 0;
    double max_time = 0;
    vector<double> part_time_list;
    vector<double> global_time_list(total_num);
    
    map<int, vector<int>> tile_box_map = calTileBox(shp_files, 0, 20);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> z_dis(0, 20);

    while (1){
        int z = z_dis(gen);
        vector<int> tile_box = tile_box_map[z];
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);
        int x = x_dis(gen);
        int y = y_dis(gen);

        
        linestringNode *tile_node = locTileNode(z, x, y, RNode);
        if (tile_node){
            double tile_time = 0;

            if (tile_node->node_level <= level_thd){
                if (tile_node->node_fcount == 0 or tile_node->node_flen == 0){
                    continue;
                }

                if (tile_node->node_fcount < count_thd){
                    if (tile_node->node_flen < len_thd){
                        tile_time = DaPGPlot(z, x, y, 2, RNode, level_thd);
                    }
                    else{
                        tile_time = DiPGPlot(z, x, y, 2, RNode);
                    }
                }
                else{
                    tile_time = DiPGPlot(z, x, y, 2, RNode);
                }
            }
            else{
                auto s_time = chrono::high_resolution_clock::now();

                dataDrivenPlot(z, x, y, 255, 0, 0, 1, 2, RNode, level_thd);

                auto e_time = chrono::high_resolution_clock::now();

                chrono::duration<double> time_use = e_time - s_time;

                tile_time = time_use.count();
            }

            part_time_list.push_back(tile_time);

            total_time += tile_time;
            tile_num++;
        }

        if (tile_num == total_num / process_num){
            break;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Allgather(part_time_list.data(), total_num / process_num, MPI_DOUBLE, 
                  global_time_list.data(), total_num / process_num, MPI_DOUBLE, 
                  MPI_COMM_WORLD);
    
    
    if (cpu_rank == 0){
        for (auto tile_time: global_time_list){
            cout << tile_time << endl;
        }

        cout << max_time << " " << total_num / max_time << endl;
    }
}


void calTimeSpeed(int cpu_rank, int process_num, vector<string> shp_files, polygonNode *RNode, short level_thd, int count_thd, double len_thd, double area_thd){
    int total_num = 5000;
    int tile_num = 0;
    double total_time = 0;
    double max_time = 0;
    vector<double> part_time_list;
    vector<double> global_time_list(total_num);

    map<int, vector<int>> tile_box_map = calTileBox(shp_files, 0, 20);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> z_dis(0, 20);

    while (1){
        int z = z_dis(gen);
        vector<int> tile_box = tile_box_map[z];
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);
        int x = x_dis(gen);
        int y = y_dis(gen);

        
        polygonNode *tile_node = locTileNode(z, x, y, RNode);
        if (tile_node && tile_node->topo_type == 'i'){
            double tile_time = 0;

            if (tile_node->node_level <= level_thd){
                if (tile_node->node_fcount == 0 or tile_node->node_flen == 0 or tile_node->node_farea == 0){
                    continue;
                }

                if (tile_node->node_fcount < count_thd){
                    if (tile_node->node_flen < len_thd){
                        if (tile_node->node_farea < area_thd){
                            tile_time = DaPGPlot(z, x, y, 2, RNode, level_thd);
                        }
                        else{
                            tile_time = DiPGPlot(z, x, y, 2, RNode);
                        }
                    }
                    else{
                        tile_time = DiPGPlot(z, x, y, 2, RNode);
                    }
                }
                else{
                    tile_time = DiPGPlot(z, x, y, 2, RNode);
                }
            }
            else{
                auto s_time = chrono::high_resolution_clock::now();

                dataDrivenPlot(z, x, y, 255, 0, 0, 1, 2, RNode, level_thd);

                auto e_time = chrono::high_resolution_clock::now();

                chrono::duration<double> time_use = e_time - s_time;

                tile_time = time_use.count();
            }

            part_time_list.push_back(tile_time);

            total_time += tile_time;
            tile_num++;
        }

        if (tile_num == total_num / process_num){
            break;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Allgather(part_time_list.data(), total_num / process_num, MPI_DOUBLE, 
                  global_time_list.data(), total_num / process_num, MPI_DOUBLE, 
                  MPI_COMM_WORLD);
    
    
    if (cpu_rank == 0){
        for (auto tile_time: global_time_list){
            cout << tile_time << endl;
        }

        cout << max_time << " " << total_num / max_time << endl;
    }
}



int main(int argc, char **argv){
    int cpu_rank, process_num; 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &cpu_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);

    string data_path = argv[1];
    string shp_type = argv[2];

    vector<string> shp_files;

    if (getShpFiles(data_path, shp_files)){
        if (shp_type == "point"){
            short level_thd = calTransiLevel(shp_files, "4326") + 9;
            int count_thd = 600;

            pointNode *pointRNode = pointIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);
        
            // calculate the total time use of PATD
            calTimeSpeed(cpu_rank, process_num, shp_files, pointRNode, level_thd, count_thd);
        }
        else if (shp_type == "linestring"){
            short level_thd = calTransiLevel(shp_files, "4326") + 9;
            int count_thd = 700;
            double len_thd = 200;

            linestringNode *linestringRNode = linestringIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the total time use of PATD
            calTimeSpeed(cpu_rank, process_num, shp_files, linestringRNode, level_thd, count_thd, len_thd);
        }
        else if (shp_type == "polygon"){
            short level_thd = calTransiLevel(shp_files, "4326") + 9;
            int count_thd = 900;
            double len_thd = 300;
            double area_thd = 3;

            polygonNode *polygonRNode = polygonIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the total time use of PATD
            calTimeSpeed(cpu_rank, process_num, shp_files, polygonRNode, level_thd, count_thd, len_thd, area_thd);
        }
    }

    MPI_Finalize();
    return 0;
}