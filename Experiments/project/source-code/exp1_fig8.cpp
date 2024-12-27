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


void calTime(int cpu_rank, int process_num, vector<string> shp_files, pointNode *pointRNode, short s_level, short e_level, vector<double> &DiPG_list, vector<double> &DaPG_list){
    for (int level = s_level; level <= e_level; level++){
        double DiPG_total_time = 0;
        double DaPG_total_time = 0;
        int tile_count = 0;

        vector<int> tile_box = calTileBox(shp_files, level);

        for (int x = tile_box[0]; x <= tile_box[2]; x++){
            for (int y = tile_box[1]; y <= tile_box[3]; y++){
                if (cpu_rank == tile_count % process_num){
                    pointNode *tile_node = locTileNode(level, x, y, pointRNode);
                    if (tile_node){
                        double DiPG_time = DiPGPlot(level, x, y, 2, pointRNode);
                        double DaPG_time = DaPGPlot(level, x, y, 2, pointRNode, e_level);

                        DiPG_total_time += DiPG_time;
                        DaPG_total_time += DaPG_time;
                    }
                }
                tile_count++;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double DiPG_max_time = 0;
        double DaPG_max_time = 0;
        MPI_Reduce(&DiPG_total_time, &DiPG_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&DaPG_total_time, &DaPG_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


        DiPG_list.push_back(DiPG_max_time);
        DaPG_list.push_back(DaPG_max_time);
    }
}


void calTime(int cpu_rank, int process_num, vector<string> shp_files, linestringNode *linestringRNode, short s_level, short e_level, vector<double> &DiPG_list, vector<double> &DaPG_list){
    for (int level = s_level; level <= e_level; level++){
        double DiPG_total_time = 0;
        double DaPG_total_time = 0;
        int tile_count = 0;

        vector<int> tile_box = calTileBox(shp_files, level);

        for (int x = tile_box[0]; x <= tile_box[2]; x++){
            for (int y = tile_box[1]; y <= tile_box[3]; y++){
                if (cpu_rank == tile_count % process_num){
                    linestringNode *tile_node = locTileNode(level, x, y, linestringRNode);
                    if (tile_node){
                        double DiPG_time = DiPGPlot(level, x, y, 2, linestringRNode);
                        double DaPG_time = DaPGPlot(level, x, y, 2, linestringRNode, e_level);

                        DiPG_total_time += DiPG_time;
                        DaPG_total_time += DaPG_time;
                    }
                }
                tile_count++;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double DiPG_max_time = 0;
        double DaPG_max_time = 0;
        MPI_Reduce(&DiPG_total_time, &DiPG_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&DaPG_total_time, &DaPG_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


        DiPG_list.push_back(DiPG_max_time);
        DaPG_list.push_back(DaPG_max_time);
    }
}


void calTime(int cpu_rank, int process_num, vector<string> shp_files, polygonNode *polygonRNode, short s_level, short e_level, vector<double> &DiPG_list, vector<double> &DaPG_list){
    for (int level = s_level; level <= e_level; level++){
        double DiPG_total_time = 0;
        double DaPG_total_time = 0;
        int tile_count = 0;

        vector<int> tile_box = calTileBox(shp_files, level);

        for (int x = tile_box[0]; x <= tile_box[2]; x++){
            for (int y = tile_box[1]; y <= tile_box[3]; y++){
                if (cpu_rank == tile_count % process_num){
                    polygonNode *tile_node = locTileNode(level, x, y, polygonRNode);
                    if (tile_node){
                        double DiPG_time = DiPGPlot(level, x, y, 3, polygonRNode);
                        double DaPG_time = DaPGPlot(level, x, y, 3, polygonRNode, e_level);

                        DiPG_total_time += DiPG_time;
                        DaPG_total_time += DaPG_time;
                    }
                }
                tile_count++;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double DiPG_max_time = 0;
        double DaPG_max_time = 0;
        MPI_Reduce(&DiPG_total_time, &DiPG_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&DaPG_total_time, &DaPG_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


        DiPG_list.push_back(DiPG_max_time);
        DaPG_list.push_back(DaPG_max_time);
    }
}


void calSpeed(int process_num, vector<string> shp_files, pointNode *pointRNode, short s_level, short e_level, vector<double> &DiPG_list, vector<double> &DaPG_list){
    int total_num = 1000;

    for (int level = s_level; level <= e_level; level++){
        int tile_num = 0;
        double DiPG_total_time = 0;
        double DaPG_total_time = 0;

        vector<int> tile_box = calTileBox(shp_files, level);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);

        while(1){
            int x = x_dis(gen);
            int y = y_dis(gen);
            
            pointNode *tile_node = locTileNode(level, x, y, pointRNode);
            if (tile_node){
                double DiPG_time = DiPGPlot(level, x, y, 2, pointRNode);
                double DaPG_time = DaPGPlot(level, x, y, 2, pointRNode, e_level);

                DiPG_total_time += DiPG_time;
                DaPG_total_time += DaPG_time;
                
                tile_num++;
            }

            if (tile_num == total_num / process_num){
                break;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        double DiPG_max_time = 0;
        double DaPG_max_time = 0;
        MPI_Reduce(&DiPG_total_time, &DiPG_max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&DaPG_total_time, &DaPG_max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        double DiPG_speed = total_num / DiPG_max_time;
        double DaPG_speed = total_num / DaPG_max_time;

        DiPG_list.push_back(DiPG_speed);
        DaPG_list.push_back(DaPG_speed);
    }
}


void calSpeed(int process_num, vector<string> shp_files, linestringNode *linestringRNode, short s_level, short e_level, vector<double> &DiPG_list, vector<double> &DaPG_list){
    int total_num = 1000;

    for (int level = s_level; level <= e_level; level++){
        int tile_num = 0;
        double DiPG_total_time = 0;
        double DaPG_total_time = 0;

        vector<int> tile_box = calTileBox(shp_files, level);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);

        while(1){
            int x = x_dis(gen);
            int y = y_dis(gen);
            
            linestringNode *tile_node = locTileNode(level, x, y, linestringRNode);
            if (tile_node){
                double DiPG_time = DiPGPlot(level, x, y, 2, linestringRNode);   
                double DaPG_time = DaPGPlot(level, x, y, 2, linestringRNode, e_level);

                DiPG_total_time += DiPG_time;
                DaPG_total_time += DaPG_time;
                
                tile_num++;
            }

            if (tile_num == total_num / process_num){
                break;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        double DiPG_max_time = 0;
        double DaPG_max_time = 0;
        MPI_Reduce(&DiPG_total_time, &DiPG_max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&DaPG_total_time, &DaPG_max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        double DiPG_speed = total_num / DiPG_max_time;
        double DaPG_speed = total_num / DaPG_max_time;

        DiPG_list.push_back(DiPG_speed);
        DaPG_list.push_back(DaPG_speed);
    }
}


void calSpeed(int process_num, vector<string> shp_files, polygonNode *polygonRNode, short s_level, short e_level, vector<double> &DiPG_list, vector<double> &DaPG_list){
    int total_num = 1000;

    for (int level = s_level; level <= e_level; level++){
        int tile_num = 0;
        double DiPG_total_time = 0;
        double DaPG_total_time = 0;

        vector<int> tile_box = calTileBox(shp_files, level);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> x_dis(tile_box[0], tile_box[2]);
        uniform_int_distribution<> y_dis(tile_box[1], tile_box[3]);

        while(1){
            int x = x_dis(gen);
            int y = y_dis(gen);
            
            polygonNode *tile_node = locTileNode(level, x, y, polygonRNode);
            if (tile_node){
                double DiPG_time = DiPGPlot(level, x, y, 3, polygonRNode);   
                double DaPG_time = DaPGPlot(level, x, y, 3, polygonRNode, e_level);

                DiPG_total_time += DiPG_time;
                DaPG_total_time += DaPG_time;
                
                tile_num++;
            }

            if (tile_num == total_num / process_num){
                break;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        double DiPG_max_time = 0;
        double DaPG_max_time = 0;
        MPI_Reduce(&DiPG_total_time, &DiPG_max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&DaPG_total_time, &DaPG_max_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        double DiPG_speed = total_num / DiPG_max_time;
        double DaPG_speed = total_num / DaPG_max_time;

        DiPG_list.push_back(DiPG_speed);
        DaPG_list.push_back(DaPG_speed);
    }
}




int main(int argc, char **argv){
    int cpu_rank, process_num; 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &cpu_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);

    string data_path = argv[1];
    string shp_type = argv[2];
    short min_level = stoi(argv[3]);
    short max_level = stoi(argv[4]);

    vector<string> shp_files;
    vector<double> DiPG_time_list, DaPG_time_list;
    vector<double> DiPG_speed_list, DaPG_speed_list;

    if (getShpFiles(data_path, shp_files)){
        if (shp_type == "point"){
            pointNode *pointRNode = pointIndex(shp_files, "4326", max_level);

            MPI_Barrier(MPI_COMM_WORLD);
        
            // calculate the time use of DiPG and DaPG in each zoom level
            calTime(cpu_rank, process_num, shp_files, pointRNode, min_level, max_level, DiPG_time_list, DaPG_time_list);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of DiPG and DaPG in each zoom level
            calSpeed(process_num, shp_files, pointRNode, min_level, max_level, DiPG_speed_list, DaPG_speed_list);
        }
        else if (shp_type == "linestring"){
            linestringNode *linestringRNode = linestringIndex(shp_files, "4326", max_level);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the time use of DiPG and DaPG in each zoom level
            calTime(cpu_rank, process_num, shp_files, linestringRNode, min_level, max_level, DiPG_time_list, DaPG_time_list);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of DiPG and DaPG in each zoom level
            calSpeed(process_num, shp_files, linestringRNode, min_level, max_level, DiPG_speed_list, DaPG_speed_list);
        }
        else if (shp_type == "polygon"){
            polygonNode *polygonRNode = polygonIndex(shp_files, "4326", max_level);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the time use of DiPG and DaPG in each zoom level
            calTime(cpu_rank, process_num, shp_files, polygonRNode, min_level, max_level, DiPG_time_list, DaPG_time_list);

            MPI_Barrier(MPI_COMM_WORLD);

            // calculate the plot speed of DiPG and DaPG in each zoom level
            calSpeed(process_num, shp_files, polygonRNode, min_level, max_level, DiPG_speed_list, DaPG_speed_list);
        }

        if (cpu_rank == 0){
            // print total time and plot speed
            for (int i = 0; i < DiPG_time_list.size(); i++){
                cout << DiPG_time_list[i] << " " << DaPG_time_list[i] << " " << DiPG_speed_list[i] << " " << DaPG_speed_list[i] << endl;
            }
        }
    }

    MPI_Finalize();
    return 0;
}