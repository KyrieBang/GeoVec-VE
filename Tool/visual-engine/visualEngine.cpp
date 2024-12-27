#include <filesystem>
#include "crow.h"
#include "png.h"
#include <shapefil.h>
#include "indexNode.hpp"
#include "buildIndex.hpp"
#include "pointRender.hpp"
#include "linestringRender.hpp"
#include "polygonRender.hpp"
#define TILE_SIZE 256
#define L	20037508.34
namespace fs = filesystem;

struct ServerMiddleware
{
	std::string message;

	ServerMiddleware()
	{
		message = "foo";
	}


	void setMessage( std::string newMsg )
	{
		message = newMsg;
	}


	struct context
	{
	};

	void before_handle( crow::request & /*req*/, crow::response & /*res*/, context & /*ctx*/ )
	{
		CROW_LOG_DEBUG << " - MESSAGE: " << message;
	}


	void after_handle( crow::request & /*req*/, crow::response & /*res*/, context & /*ctx*/ )
	{
		/* no-op */
	}
};


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


bool checkShpFiles(vector<string> shp_files, string &data_espg){
    string espg_tmp;
    
    for (int i = 0; i < shp_files.size(); i++){
        string shp_path = shp_files[i];

        // check the srs of each shapefile
        string prj_path = shp_path.substr(0, shp_path.find(".shp")) + ".prj";
        ifstream prjFile(prj_path);
        if (prjFile.is_open()) {
            string prjContent((istreambuf_iterator<char>(prjFile)), istreambuf_iterator<char>());

            if (prjContent.find("PROJCS") == string::npos and prjContent.find("WGS_1984") != string::npos){
                data_espg = "4326";
            }
            else if (prjContent.find("PROJCS") != string::npos and prjContent.find("Web_Mercator") != string::npos){
                data_espg = "3857";
            }
            else{
                cout << "【ERROR】 please input the shapefile with ESPG is 4326 or 3857 !" << endl;
                cout << "【ERROR】 " << prj_path << endl;
                return false;
            }

            if (i == 0){
                espg_tmp = data_espg;
            }
            
            if (espg_tmp != data_espg){
                cout << "【ERROR】 please check wheather the ESPG of the shapefiles are same !" << endl;
                for (int j = 0; j <= i; j++){
                    cout << "【ERROR】 " << shp_files[j].substr(0, shp_files[j].find(".shp")) + ".prj" << endl;
                }
                return false;
            }
        }
        else{
            cout << "【ERROR】 please check the .prj file whether exists !" << endl;
            cout << "【ERROR】 " << prj_path << endl;
            return false;
        }
    }

    return true; 
}




int main(int argc, char **argv){    
    // get parameters
    string point_data_path = argv[1];
    string linestring_data_path = argv[2];
    string polygon_data_path = argv[3];
    vector<string> point_shp_files, linestring_shp_files, polygon_shp_files;
    string point_espg, linestring_espg, polygon_espg;
    short point_transiLevel, linestring_transiLevel, polygon_transiLevel;
    pointNode *pointRNode;
    linestringNode *linestringRNode;
    polygonNode *polygonRNode;


    // ********************************************************************************************************************* //
	// read point shapefile
	// ********************************************************************************************************************* //
    if (point_data_path == "null"){
        pointRNode = nullptr;
    }
    else{
        if (!getShpFiles(point_data_path, point_shp_files) or !checkShpFiles(point_shp_files, point_espg)){
            return 0;
        }

        cout << "--------------------------------------------------------------------------------" << endl;
        point_transiLevel = calTransiLevel(point_shp_files, point_espg) + 9;
        pointRNode = pointIndex(point_shp_files, point_espg, point_transiLevel);
        cout << "--------------------------------------------------------------------------------" << endl;
    }

    // ********************************************************************************************************************* //
	// read linestring shapefile
	// ********************************************************************************************************************* //
    if (linestring_data_path == "null"){
        linestringRNode = nullptr;
    }
    else{
        if (!getShpFiles(linestring_data_path, linestring_shp_files) or !checkShpFiles(linestring_shp_files, linestring_espg)){
            return 0;
        }

        cout << "--------------------------------------------------------------------------------" << endl;
        linestring_transiLevel = calTransiLevel(linestring_shp_files, linestring_espg) + 9;
        linestringRNode = linestringIndex(linestring_shp_files, linestring_espg, linestring_transiLevel);
        cout << "--------------------------------------------------------------------------------" << endl;
    }

    // ********************************************************************************************************************* //
	// read polygon shapefile
	// ********************************************************************************************************************* //
    if (polygon_data_path == "null"){
        polygonRNode = nullptr;
    }
    else{
        if (!getShpFiles(polygon_data_path, polygon_shp_files) or !checkShpFiles(polygon_shp_files, polygon_espg)){
            return 0;
        }

        cout << "--------------------------------------------------------------------------------" << endl;
        polygon_transiLevel = calTransiLevel(polygon_shp_files, polygon_espg) + 9;
        polygonRNode = polygonIndex(polygon_shp_files, polygon_espg, polygon_transiLevel);
        cout << "--------------------------------------------------------------------------------" << endl;
    }

    cout << "<<<<<<< read shapefile successfully ! >>>>>>>" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    

    // ********************************************************************************************************************* //
	// WTMS of vector data visualization 
	// ********************************************************************************************************************* //
    crow::App<ServerMiddleware> app;
    CROW_ROUTE(app, "/GeoVec-VE/<string>/<int>/<int>/<int>/<float>/<float>/<int>/<int>/<int>.png").name("GeoVec-VE")
		([&] (const crow::request & req, crow::response & res, string shpType, int R, int G, int B, float AD, float width, int z, int x, int y) {
			string img_str;
            try{
                if (shpType == "point"){
                    if (pointRNode){
                        if (z <= point_transiLevel){
                            img_str = disDrivenPlot(z, x, y, R, G, B, AD, width, pointRNode);
                        }
                        else{
                            img_str = dataDrivenPlot(z, x, y, R, G, B, AD, width, pointRNode, point_transiLevel);
                        }
                    }
                }
                else if (shpType == "line"){
                    if (linestringRNode){
                        if (z <= linestring_transiLevel){
                            img_str = disDrivenPlot(z, x, y, R, G, B, AD, width, linestringRNode);
                        }
                        else{
                            img_str = dataDrivenPlot(z, x, y, R, G, B, AD, width, linestringRNode, linestring_transiLevel);
                        }
                    }
                }
                else if (shpType == "polygon"){
                    if (polygonRNode){
                        if (z <= polygon_transiLevel){
                            img_str = disDrivenPlot(z, x, y, R, G, B, AD, width, polygonRNode);
                        }
                        else{
                            img_str = dataDrivenPlot(z, x, y, R, G, B, AD, width, polygonRNode, polygon_transiLevel);
                        }
                    }
                }

				res.write(img_str);
				res.set_header("Content-Type", "image/png" );
				res.set_header("Access-Control-Allow-Origin", "*");
				res.end();
            }
            catch(const char* msg){
				res.code = 500;
                ostringstream os;
				os << msg << "\n";
				res.write(os.str());
				res.set_header("Content-Type", "text/html");
				res.end();
            }
        });
    app.port(10085).multithreaded().run();
}