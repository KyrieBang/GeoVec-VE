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


void calTime(int cpu_rank, int process_num, vector<string> shp_files, pointNode *RNode, short level_thd, int count_thd){
    double total_time = 0;
    double max_time = 0;
    pointNode *tile_node;
    
    for (int z = 0; z <= 8; z++){
        int tile_count = 0;

        vector<int> tile_box = calTileBox(shp_files, z);

        for (int x = tile_box[0]; x <= tile_box[2]; x++){
            for (int y = tile_box[1]; y <= tile_box[3]; y++){
                if (cpu_rank == tile_count % process_num){
                    tile_node = locTileNode(z, x, y, RNode);
                    if (tile_node){
                        total_time += DiPGPlot(z, x, y, 2, RNode);
                    }
                }
                tile_count++;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (cpu_rank == 0){
        cout << max_time << endl;
    }
}


void calTime(int cpu_rank, int process_num, vector<string> shp_files, linestringNode *RNode, short level_thd, int count_thd, double len_thd){
    double total_time = 0;
    double max_time = 0;
    linestringNode *tile_node;

    for (int z = 0; z <= 8; z++){
        int tile_count = 0;

        vector<int> tile_box = calTileBox(shp_files, z);

        for (int x = tile_box[0]; x <= tile_box[2]; x++){
            for (int y = tile_box[1]; y <= tile_box[3]; y++){
                if (cpu_rank == tile_count % process_num){
                    tile_node = locTileNode(z, x, y, RNode);
                    if (tile_node){
                        total_time += DiPGPlot(z, x, y, 2, RNode);
                    }
                }
                tile_count++;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (cpu_rank == 0){
        cout << max_time << endl;
    }
}


void calTime(int cpu_rank, int process_num, vector<string> shp_files, polygonNode *RNode, short level_thd, int count_thd, double len_thd, double area_thd){
    double total_time = 0;
    double max_time = 0;
    polygonNode *tile_node;
    
    for (int z = 0; z <= 8; z++){
        int tile_count = 0;

        vector<int> tile_box = calTileBox(shp_files, z);

        for (int x = tile_box[0]; x <= tile_box[2]; x++){
            for (int y = tile_box[1]; y <= tile_box[3]; y++){
                if (cpu_rank == tile_count % process_num){
                    tile_node = locTileNode(z, x, y, RNode);
                    if (tile_node){
                        total_time += DiPGPlot(z, x, y, 2, RNode);
                    }
                }
                tile_count++;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (cpu_rank == 0){
        cout << max_time << endl;
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
            short level_thd = calTransiLevel(shp_files, "4326") + 8;
            pointNode *pointRNode = pointIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);
        
            // calculate the total time in zoom levels 0-8
            calTime(cpu_rank, process_num, shp_files, pointRNode, level_thd, 600);
        }
        else if (shp_type == "linestring"){
            short level_thd = calTransiLevel(shp_files, "4326") + 8;
            linestringNode *linestringRNode = linestringIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the total time in zoom levels 0-8
            calTime(cpu_rank, process_num, shp_files, linestringRNode, level_thd, 700, 200);
        }
        else if (shp_type == "polygon"){
            short level_thd = calTransiLevel(shp_files, "4326") + 8;
            polygonNode *polygonRNode = polygonIndex(shp_files, "4326", level_thd);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the total time in zoom levels 0-8
            calTime(cpu_rank, process_num, shp_files, polygonRNode, level_thd, 900, 300, 3);
        }
    }

    MPI_Finalize();
    return 0;
}