#include "png.h"
#include "linestringRender.hpp"
#include <omp.h>
#include <mpi.h>
#include <chrono>
#include <mapnik/map.hpp>
#include <mapnik/agg_renderer.hpp>
#include <mapnik/image.hpp>
#include <mapnik/image_util.hpp>
#include <mapnik/layer.hpp>
#include <mapnik/rule.hpp>
#include <mapnik/feature_type_style.hpp>
#include <mapnik/symbolizer.hpp>
#include <mapnik/text/placements/dummy.hpp>
#include <mapnik/text/text_properties.hpp>
#include <mapnik/text/formatting/text.hpp>
#include <mapnik/datasource_cache.hpp>
#include <mapnik/datasource.hpp>
#include <mapnik/font_engine_freetype.hpp>
#include <mapnik/expression.hpp>
#include <mapnik/color_factory.hpp>
#include <mapnik/unicode.hpp>
#include <mapnik/save_map.hpp>
#include <mapnik/cairo_io.hpp>
#include <mapnik/geometry.hpp>
#include <mapnik/geometry_envelope.hpp>
#include <mapnik/query.hpp>
#include <mapnik/geometry_centroid.hpp>
using namespace mapnik;
using namespace mapnik::geometry;
#define L 20037508.3427892
#define TILE_SIZE 256


const bool mapnik_statue = datasource_cache::instance().register_datasources("/usr/lib/mapnik/3.0/input");
const string srs_merc = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0.0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs +over";


// determine whether the node correspongding to the pixel exists
bool nodeIfExist(int z, float tile_box[], float pix_x, float pix_y, linestringNode *seek){
    float min_x = tile_box[0];
    float min_y = tile_box[1];
    float max_x = tile_box[2];
    float max_y = tile_box[3];
    float mid_x, mid_y;
    int divided_num = 0;

    while (divided_num < 8){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // pixel in LU sub-region
        if (pix_x < mid_x && pix_y >= mid_y){
            if (seek->LUNode == nullptr){
                return false;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // pixel in RU sub-region
        else if (pix_x >= mid_x && pix_y >= mid_y){
            if (seek->RUNode == nullptr){
                return false;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // pixel in LB sub-region
        else if (pix_x < mid_x && pix_y < mid_y){
            if (seek->LBNode == nullptr){
                return false;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // pixel in RB sub-region
        else if (pix_x >= mid_x && pix_y < mid_y){
            if (seek->RBNode == nullptr){
                return false;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }else {
            return false;
        }
        divided_num ++;
    }
    return true;
}


// determine whether the node correspongding to the tile<z,x,y> exists
linestringNode * locTileNode(int z, int x, int y, linestringNode *seek){
    float tile_Rz = L / pow(2, z-1);
    float tile_x = (x + 0.5) * tile_Rz - L;
    float tile_y = L - (y + 0.5) * tile_Rz;
    float min_x = -L;
    float max_x = L;
    float min_y = -L;
    float max_y = L;
    float mid_x, mid_y;
    int divided_num = 0;

    while (divided_num < z){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // tile_center_pt in LU sub-region
        if (tile_x < mid_x && tile_y >= mid_y){
            if (seek->LUNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in RU sub-region
        else if (tile_x >= mid_x && tile_y >= mid_y){
            if (seek->RUNode == nullptr){
                seek = nullptr;
                break;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in LB sub-region
        else if (tile_x < mid_x && tile_y < mid_y){
            if (seek->LBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // tile_center_pt in RB sub-region
        else if (tile_x >= mid_x && tile_y < mid_y){
            if (seek->RBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }
        divided_num ++;
    }
    return seek;
}


// find the node which links to rtree
linestringNode * locRtreeNode(int did_num, int tile_x, int tile_y, linestringNode *seek, short switch_level){
    float min_x = -L;
    float max_x = L;
    float min_y = -L;
    float max_y = L;
    float mid_x, mid_y;
    int divided_num = 0;
    linestringNode *middle_node = nullptr;

    while (divided_num < did_num){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // tile_center_pt in LU sub-region
        if (tile_x < mid_x && tile_y >= mid_y){
            if (seek->LUNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in RU sub-region
        else if (tile_x >= mid_x && tile_y >= mid_y){
            if (seek->RUNode == nullptr){
                seek = nullptr;
                break;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in LB sub-region
        else if (tile_x < mid_x && tile_y < mid_y){
            if (seek->LBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // tile_center_pt in RB sub-region
        else if (tile_x >= mid_x && tile_y < mid_y){
            if (seek->RBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }
        
        // get the node which links to r-tree
        if (divided_num == switch_level){
            middle_node = seek;
        }

        divided_num ++;
    }

    if (seek){
        return middle_node;
    }
    return seek;
}


// determine whether the node correspongding to the nearest pixel exists
linestringNode * locPixNode(int z, float tile_box[], float pix_x, float pix_y, linestringNode *seek, linestringNode *RNode){
    float min_x = tile_box[0];
    float min_y = tile_box[1];
    float max_x = tile_box[2];
    float max_y = tile_box[3];
    float mid_x, mid_y;

    if (abs(pix_x) <= L && abs(pix_y) <= L){
        if (pix_x < min_x || pix_x > max_x || pix_y < min_y || pix_y > max_y){
            seek = RNode;
            min_x = -L;
            max_x = L;
            min_y = -L;
            max_y = L;
        }
    }
    else{
        seek = nullptr;
        return seek;
    }

    int level = 8 + z - seek->node_level;
    int divided_num = 0;

    while (divided_num < level){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // pixel in LU sub-region
        if (pix_x < mid_x && pix_y >= mid_y){
            if (seek->LUNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // pixel in RU sub-region
        else if (pix_x >= mid_x && pix_y >= mid_y){
            if (seek->RUNode == nullptr){
                seek = nullptr;
                break;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // pixel in LB sub-region
        else if (pix_x < mid_x && pix_y < mid_y){
            if (seek->LBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // pixel in RB sub-region
        else if (pix_x >= mid_x && pix_y < mid_y){
            if (seek->RBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }
        divided_num++;
    }
    return seek;
}


// display-driven tile plotting
void dis_tileRender(int z, int x, int y, int R, int G, int B, float AD, float width, linestringNode *tile_node, linestringNode *linestringRNode, png_bytep *png_ptr){
    float Rz = L / (128 << z);
    float R0 = width * Rz; 
    float R1 = R0 - sqrt(2) * Rz / 4.0;
    float R2 = R0 + sqrt(2) * Rz / 4.0;
    float Rz4 = 0.25 * Rz;
    int level = int(width + sqrt(2) / 4.0 - 0.5);
    // calcalate the BBox of tile <z,x,y>
    float tile_Rz = L / pow(2, z-1);
    float tile_xMin = x * tile_Rz - L;
    float tile_xMax = (x + 1) * tile_Rz - L;
    float tile_yMin = L - (y + 1) * tile_Rz;
    float tile_yMax = L - y * tile_Rz;
    float tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            int tile_index = i * 256 + j;
            short final_value = 0;
            // calculate center coords of each pixel
            float pix_x = (256 * x + j + 0.5) * Rz - L;
            float pix_y = L - (256 * y + i + 0.5) * Rz;
            // determine whether the node correspongding to the pixel exists
            if (nodeIfExist(z, tile_box, pix_x, pix_y, tile_node)){
                png_ptr[i][4 * j] = R;
                png_ptr[i][4 * j + 1] = G;
                png_ptr[i][4 * j + 2] = B;
                png_ptr[i][4 * j + 3] = AD * final_value * 64 - 1;
                continue;
            }
            short pix_value = 0;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (m == -l || m == l || n == -l || n == l){
                            linestringNode *node_floor = locPixNode(z, tile_box, pix_x + m * Rz, pix_y + n * Rz, tile_node, linestringRNode);
                            if (node_floor != nullptr){
                                if (Rz * sqrt(m * m + n * n) < R1){
                                    final_value = 4;
                                    break;
                                }
                                else if (Rz * sqrt(m * m + n * n) < R2){ 
                                    if (sqrt(pow(Rz4 - m * Rz, 2) + pow(Rz4 - n * Rz, 2)) < R0){
                                        pix_value += 1;
                                    }
                                    if (sqrt(pow(- Rz4 - m * Rz, 2) + pow(- Rz4 - n * Rz, 2)) < R0){
                                        pix_value += 1;
                                    }
                                    if (sqrt(pow(Rz4 - m * Rz, 2) + pow(- Rz4 - n * Rz, 2)) < R0){
                                        pix_value += 1;
                                    }  
                                    if (sqrt(pow(- Rz4 - m * Rz, 2) + pow(Rz4 - n * Rz, 2)) < R0){
                                        pix_value += 1; 
                                    }

                                    if (pix_value == 4){
                                        final_value = 4;
                                        break;
                                    }
                                    else if (pix_value > final_value){
                                        final_value = pix_value;
                                    }
                                    pix_value = 0; 
                                }
                            }
                        }
                        
                    }
                    if (final_value == 4)
                        break;
                }
                if (final_value == 4)
                    break;
            }
            // generate pixel value
            if (final_value > 0){
                png_ptr[i][4 * j] = R;
                png_ptr[i][4 * j + 1] = G;
                png_ptr[i][4 * j + 2] = B;
                png_ptr[i][4 * j + 3] = AD * final_value * 64 - 1;
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
        }
    }
}


// data-driven tile plotting
string data_tileRender_l(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, int R, int G, int B, float AD, float width, datasource_ptr line_ds){
    Map map(256, 256);
    map.set_srs(srs_merc);

    feature_type_style line_style;
    {
        rule r;
        {
            line_symbolizer line_sym;
            put(line_sym, keys::stroke, color(R, G, B));
            put(line_sym, keys::stroke_width, width);
            put(line_sym, keys::stroke_opacity, AD);
            put(line_sym, keys::smooth, 1.0);
            r.append(std::move(line_sym));
        }
        line_style.add_rule(std::move(r));
    }
    map.insert_style("line_style", std::move(line_style));

    {
        layer lyr("line");
        lyr.set_datasource(line_ds);
        lyr.add_style("line_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    map.zoom_to_box(box2d<double>(tile_xMin, tile_yMin, tile_xMax, tile_yMax));
    image_rgba8 im(256, 256);
    agg_renderer<image_rgba8> ren(map, im);
    ren.apply();
    ostringstream im_os;
    save_to_stream(im, im_os, "png");
    return im_os.str();
}


// diaplay-driven visualization
string disDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, linestringNode *linestringRNode){
    png_bytep * tile_ptr = (png_bytep *) malloc(256 * sizeof(png_bytep));
    for (int n = 0; n < TILE_SIZE; n++)
        tile_ptr[n] = (png_bytep) malloc(1024);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr); 
    FILE *temp_png = tmpfile();
    png_init_io(png_ptr, temp_png); 
    png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    linestringNode *tile_node = locTileNode(z, x, y, linestringRNode);
    if (tile_node != nullptr){
        dis_tileRender(z, x, y, R, G, B, AD, width, tile_node, linestringRNode, tile_ptr);
    }
    else{
        # pragma omp parallel for num_threads(4)
        for (int i = 0; i < 256; i++){
            for (int j = 0; j < 256 * 4; j++){
                tile_ptr[i][j] = 0; 
            }
        }
    }

    png_write_image(png_ptr, tile_ptr);   
    png_write_end(png_ptr, NULL);
    for (int k = 0; k < 256; k++)
        free(tile_ptr[k]);
    free(tile_ptr);
    vector<char> pos;
    fseek(temp_png, 0, SEEK_END);
    long size = ftell(temp_png);
    rewind(temp_png);
    pos.resize(size);
    fread(&pos[0], 1, size, temp_png);
    fclose(temp_png);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    string pos_tostr = string(pos.begin(), pos.end());
    return pos_tostr;
}


// data-driven visualization
string dataDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, linestringNode *linestringRNode, short switch_level){
    string img_str;
    linestringNode *rtree_node = nullptr;
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    float tile_xCenter = tile_xMin + (tile_xMax - tile_xMin) / 2;
    float tile_yCenter = tile_yMin + (tile_yMax - tile_yMin) / 2;
    // determine whether the tile needs to be drawn
    short linestring_treeLevel = switch_level + 8;
    if (z <= linestring_treeLevel){
        rtree_node = locRtreeNode(z, tile_xCenter, tile_yCenter, linestringRNode, switch_level);
    }
    else{
        rtree_node = locRtreeNode(linestring_treeLevel, tile_xCenter, tile_yCenter, linestringRNode, switch_level);
    }
    // non-blank tile
    if (rtree_node != nullptr){
        // get the features by spatial range retrieve
        vector<BOX_PAIR> lines;
        Box box(POINT(tile_xMin, tile_yMin), POINT(tile_xMax, tile_yMax));
        rtree_node->Rtree.query(bgi::intersects(box), back_inserter(lines));
        string data_string = "";
        for (int i = 0; i < lines.size(); i++){
            XY_VEC line_cord = lines[i].second;
            data_string += "\"LINESTRING(";
            for (int j = 0; j < line_cord.size()-2; j+=2){
                data_string += (to_string(line_cord[j]) + " " + to_string(line_cord[j+1]) + ",");
            }
            data_string += (to_string(line_cord[line_cord.size()-2]) + " " + to_string(line_cord.back()) + ")\"\n");
        }
        // plot features
        parameters params;
        params["type"] = "csv";
        params["inline"] = "wkt\n" + data_string;
        datasource_ptr line_ds = datasource_cache::instance().create(params);
        img_str = data_tileRender_l(tile_xMin, tile_yMin, tile_xMax, tile_yMax, R, G, B, AD, width, line_ds);
    }
    // blank tile
    else{
        Map map(256, 256);
        const color back_color(255, 255, 255, 0);
        map.set_background(back_color);
        image_rgba8 im(256, 256);
        agg_renderer<image_rgba8> ren(map, im);
        ren.apply();
        ostringstream im_os;
        save_to_stream(im, im_os, "png");
        img_str = im_os.str();
    }
    return img_str;
}


double DiPGPlot(int z, int x, int y, float width, linestringNode *RNode){
    auto s_time = chrono::high_resolution_clock::now();
    
    png_bytep * tile_ptr = (png_bytep *) malloc(256 * sizeof(png_bytep));
    for (int n = 0; n < TILE_SIZE; n++)
        tile_ptr[n] = (png_bytep) malloc(1024);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr); 
    FILE *temp_png = tmpfile();
    png_init_io(png_ptr, temp_png); 
    png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    linestringNode *tile_node = locTileNode(z, x, y, RNode);
    dis_tileRender(z, x, y, 255, 0, 0, 1, width, tile_node, RNode, tile_ptr);

    png_write_image(png_ptr, tile_ptr);   
    png_write_end(png_ptr, NULL);

    for (int k = 0; k < 256; k++)
        free(tile_ptr[k]);
    free(tile_ptr);

    fclose(temp_png);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    auto e_time = chrono::high_resolution_clock::now();

    chrono::duration<double> time_use = e_time - s_time;

    return time_use.count();
}


double DaPGPlot(int z, int x, int y, float width, linestringNode *RNode, short transi_level){
    string data_string = "";
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;

    linestringNode *tile_node = locTileNode(z, x, y, RNode);

    vector<linestringRtree> Rtree_list;
    tile_node->locRtree(Rtree_list, transi_level);

    Box box(POINT(tile_xMin, tile_yMin), POINT(tile_xMax, tile_yMax));

    for (auto rtree: Rtree_list){
        vector<BOX_PAIR> lines;
        rtree.query(bgi::within(box), back_inserter(lines));

        for (int i = 0; i < lines.size(); i++){
            XY_VEC line_cord = lines[i].second;
            data_string += "\"LINESTRING(";
            for (int j = 0; j < line_cord.size()-2; j+=2){
                data_string += (to_string(line_cord[j]) + " " + to_string(line_cord[j+1]) + ",");
            }
            data_string += (to_string(line_cord[line_cord.size()-2]) + " " + to_string(line_cord.back()) + ")\"\n");
        }
    }

    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt\n" + data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(params);

    auto s_time = chrono::high_resolution_clock::now();

    // plot features
    string img_str = data_tileRender_l(tile_xMin, tile_yMin, tile_xMax, tile_yMax, 255, 0, 0, 1, width, line_ds);

    auto e_time = chrono::high_resolution_clock::now();

    chrono::duration<double> time_use = e_time - s_time;

    return time_use.count();
}