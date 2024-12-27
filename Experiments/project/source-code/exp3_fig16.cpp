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


vector<int> calTileBox(vector<string> shp_files, short level){
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

    double tile_Rz = L / pow(2, level-1);
    int tile_xMin = floor((data_xmin + L) / tile_Rz);
    int tile_xMax = floor((data_xmax + L) / tile_Rz);
    int tile_yMin = floor((L - data_ymax) / tile_Rz);
    int tile_yMax = floor((L - data_ymin) / tile_Rz);
    

    return {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
}


void calSpeed(int cpu_rank, int process_num, vector<string> shp_files, pointNode *RNode, short level_thd, int count_thd){
    int total_num = 1000;

    for (int z = 1; z <= 9; z+=2){
        int tile_num = 0;
        double total_time = 0;
        double max_time = 0;

        vector<int> tile_box = calTileBox(shp_files, z);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);

        while(1){
            int x = x_dis(gen);
            int y = y_dis(gen);
            
            pointNode *tile_node = locTileNode(z, x, y, RNode);
            if (tile_node){
                if (tile_node->node_fcount < count_thd){
                    total_time += DaPGPlot(z, x, y, 2, RNode, level_thd);
                }
                else{
                    total_time += DiPGPlot(z, x, y, 2, RNode);
                }

                tile_num++;
            }

            if (tile_num == total_num / process_num){
                break;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        if (cpu_rank == 0){
            cout << total_num / max_time << endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}


void calSpeed(int cpu_rank, int process_num, vector<string> shp_files, linestringNode *RNode, short level_thd, int count_thd, double len_thd){
    int total_num = 1000;

    for (int z = 1; z <= 9; z+=2){
        int tile_num = 0;
        double total_time = 0;
        double max_time = 0;

        vector<int> tile_box = calTileBox(shp_files, z);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);

        while(1){
            int x = x_dis(gen);
            int y = y_dis(gen);
            
            linestringNode *tile_node = locTileNode(z, x, y, RNode);
            if (tile_node){
                if (tile_node->node_fcount < count_thd){
                    if (tile_node->node_flen < len_thd){
                        total_time += DaPGPlot(z, x, y, 2, RNode, level_thd);
                    }
                    else{
                        total_time += DiPGPlot(z, x, y, 2, RNode);
                    }
                }
                else{
                    total_time += DiPGPlot(z, x, y, 2, RNode);
                }
                
                tile_num++;
            }

            if (tile_num == total_num / process_num){
                break;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        if (cpu_rank == 0){
            cout << total_num / max_time << endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}


void calSpeed(int cpu_rank, int process_num, vector<string> shp_files, polygonNode *RNode, short level_thd, int count_thd, double len_thd, double area_thd){
    int total_num = 1000;

    for (int z = 1; z <= 9; z+=2){
        int tile_num = 0;
        double total_time = 0;
        double max_time = 0;

        vector<int> tile_box = calTileBox(shp_files, z);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);

        while(1){
            int x = x_dis(gen);
            int y = y_dis(gen);
            
            polygonNode *tile_node = locTileNode(z, x, y, RNode);
            if (tile_node){
                if (tile_node->node_fcount < count_thd){
                    if (tile_node->node_flen < len_thd){
                        if (tile_node->node_farea < area_thd){
                            total_time += DaPGPlot(z, x, y, 2, RNode, level_thd);
                        }
                        else{
                            total_time += DiPGPlot(z, x, y, 2, RNode);
                        }
                    }
                    else{
                        total_time += DiPGPlot(z, x, y, 2, RNode);
                    }
                }
                else{
                    total_time += DiPGPlot(z, x, y, 2, RNode);
                }
                
                tile_num++;
            }

            if (tile_num == total_num / process_num){
                break;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        if (cpu_rank == 0){
            cout << total_num / max_time << endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
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
            pointNode *pointRNode = pointIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of PATD in zoom levels 1,3,5,7,9
            calSpeed(cpu_rank, process_num, shp_files, pointRNode, level_thd, 600);
        }
        else if (shp_type == "linestring"){
            short level_thd = calTransiLevel(shp_files, "4326") + 9;
            linestringNode *linestringRNode = linestringIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of PATD in zoom levels 1,3,5,7,9
            calSpeed(cpu_rank, process_num, shp_files, linestringRNode, level_thd, 700, 200);
        }
        else if (shp_type == "polygon"){
            short level_thd = calTransiLevel(shp_files, "4326") + 9;
            polygonNode *polygonRNode = polygonIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of PATD in zoom levels 1,3,5,7,9
            calSpeed(cpu_rank, process_num, shp_files, polygonRNode, level_thd, 900, 300, 3);
        }
    }

    MPI_Finalize();
    return 0;
}