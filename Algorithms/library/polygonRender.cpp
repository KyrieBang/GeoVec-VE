#include "png.h"
#include "polygonRender.hpp"
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
char nodeIfExist(int z, float tile_box[], float pix_x, float pix_y, polygonNode *seek){
    float min_x = tile_box[0];
    float min_y = tile_box[1];
    float max_x = tile_box[2];
    float max_y = tile_box[3];
    double mid_x, mid_y;

    int divided_num = 0;

    while (divided_num < 8){
        if (seek->topo_type == 'w'){
            return 'w';
        }
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // pixel in LU sub-region
        if (pix_x < mid_x && pix_y >= mid_y){
            if (seek->LUNode == nullptr){
                return 'n';
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // pixel in RU sub-region
        else if (pix_x >= mid_x && pix_y >= mid_y){
            if (seek->RUNode == nullptr){
                return 'n';
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // pixel in LB sub-region
        else if (pix_x < mid_x && pix_y < mid_y){
            if (seek->LBNode == nullptr){
                return 'n';
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // pixel in RB sub-region
        else if (pix_x >= mid_x && pix_y < mid_y){
            if (seek->RBNode == nullptr){
                return 'n';
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }
        divided_num++;
    }
    
    return seek->topo_type;
}


// determine whether the node correspongding to the tile<z,x,y> exists
polygonNode * locTileNode(int z, int x, int y, polygonNode *seek){
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
        if (seek->topo_type == 'w'){
            return seek;
        }

        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // tile_center_pt in LU sub-region
        if (tile_x < mid_x && tile_y >= mid_y){
            if (seek->LUNode == nullptr){
                return nullptr;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in RU sub-region
        else if (tile_x >= mid_x && tile_y >= mid_y){
            if (seek->RUNode == nullptr){
                return nullptr;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in LB sub-region
        else if (tile_x < mid_x && tile_y < mid_y){
            if (seek->LBNode == nullptr){
                return nullptr;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // tile_center_pt in RB sub-region
        else if (tile_x >= mid_x && tile_y < mid_y){
            if (seek->RBNode == nullptr){
                return nullptr;
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
polygonNode * locRtreeNode(int did_num, int tile_x, int tile_y, polygonNode *seek, short switch_level){
    float min_x = -L;
    float max_x = L;
    float min_y = -L;
    float max_y = L;
    float mid_x, mid_y;
    int divided_num = 0;
    polygonNode *middle_node = nullptr;

    while (divided_num < did_num){
        if (seek->topo_type == 'w'){
            return seek;
        }

        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // tile_center_pt in LU sub-region
        if (tile_x < mid_x && tile_y >= mid_y){
            if (seek->LUNode == nullptr){
                return nullptr;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in RU sub-region
        else if (tile_x >= mid_x && tile_y >= mid_y){
            if (seek->RUNode == nullptr){
                return nullptr;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in LB sub-region
        else if (tile_x < mid_x && tile_y < mid_y){
            if (seek->LBNode == nullptr){
                return nullptr;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // tile_center_pt in RB sub-region
        else if (tile_x >= mid_x && tile_y < mid_y){
            if (seek->RBNode == nullptr){
                return nullptr;
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

    if (seek->topo_type == 'w'){
        return seek;
    }

    if (seek){
        return middle_node;
    }
    return seek;
}


// determine whether the node correspongding to the nearest pixel exists
polygonNode * locPixNode(int z, float tile_box[], float pix_x, float pix_y, polygonNode *seek, polygonNode *rootNode){
    float min_x = tile_box[0];
    float min_y = tile_box[1];
    float max_x = tile_box[2];
    float max_y = tile_box[3];
    float mid_x, mid_y;

    if (abs(pix_x) <= L && abs(pix_y) <= L){
        if (pix_x < min_x || pix_x > max_x || pix_y < min_y || pix_y > max_y){
            seek = rootNode;
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
                return nullptr;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // pixel in RU sub-region
        else if (pix_x >= mid_x && pix_y >= mid_y){
            if (seek->RUNode == nullptr){
                return nullptr;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // pixel in LB sub-region
        else if (pix_x < mid_x && pix_y < mid_y){
            if (seek->LBNode == nullptr){
                return nullptr;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // pixel in RB sub-region
        else if (pix_x >= mid_x && pix_y < mid_y){
            if (seek->RBNode == nullptr){
                return nullptr;
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
void dis_tileRender(int z, int x, int y, int R, int G, int B, float AD, float width, polygonNode *tile_node, polygonNode *polygonRNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    double R0 = width * Rz; 
    double R1 = R0 - sqrt(2) * Rz / 4.0;
    double R2 = R0 + sqrt(2) * Rz / 4.0;
    double Rz4 = 0.25 * Rz;
    short level = int(width + sqrt(2) / 4.0 - 0.5);
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
            short final_value = 0;
            // calculate center coords of each pixel
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            // determine whether the node correspongding to the pixel exists
            char node_topoType = nodeIfExist(z, tile_box, pix_x, pix_y, tile_node);
            if (node_topoType != 'n'){
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
                        if (abs(m) == l || abs(n) == l){
                            double nearPix_x = pix_x + m * Rz;
                            double nearPix_y = pix_y + n * Rz;
                            polygonNode *node_floor = locPixNode(z, tile_box, nearPix_x, nearPix_y, tile_node, polygonRNode);
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
string data_tileRender_a(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, int R, int G, int B, float AD, float width, datasource_ptr polygon_ds){
    Map map(256, 256);
    map.set_srs(srs_merc);
    
    feature_type_style polygon_style;
    {
        rule r;
        {
            polygon_symbolizer poly_sym;
            put(poly_sym, keys::fill, color(R, G, B));
            put(poly_sym, keys::opacity, AD);
            put(poly_sym, keys::fill_opacity, AD);
            put(poly_sym, keys::gamma, 1.0);
            r.append(std::move(poly_sym));
        }
        polygon_style.add_rule(std::move(r));
    }
    map.insert_style("polygon_style", std::move(polygon_style));

    {
        layer lyr("polygon");
        lyr.set_datasource(polygon_ds);
        lyr.add_style("polygon_style");
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
string disDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, polygonNode *polygonRNode){
    png_bytep * tile_ptr = (png_bytep *) malloc(256 * sizeof(png_bytep));
    for (int n = 0; n < TILE_SIZE; n++)
        tile_ptr[n] = (png_bytep) malloc(1024);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr); 
    FILE *temp_png = tmpfile();
    png_init_io(png_ptr, temp_png); 
    png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    polygonNode *tile_node = locTileNode(z, x, y, polygonRNode);
    if (tile_node != nullptr){
        if (tile_node->topo_type == 'i'){
            dis_tileRender(z, x, y, R, G, B, AD, width, tile_node, polygonRNode, tile_ptr);
        }
        else{
            # pragma omp parallel for num_threads(4)
            for (int i = 0; i < 256; i++){
                for (int j = 0; j < 256; j++){
                    tile_ptr[i][4 * j] = R;
                    tile_ptr[i][4 * j + 1] = G;
                    tile_ptr[i][4 * j + 2] = B;
                    tile_ptr[i][4 * j + 3] = AD * 256 - 1;
                }
            }
        }
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
string dataDrivenPlot(int z , int x, int y, int R, int G, int B, float AD, float width, polygonNode *polygonRNode, short switch_level){
    string img_str;
    polygonNode *rtree_node = nullptr;
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    float tile_xCenter = tile_xMin + (tile_xMax - tile_xMin) / 2;
    float tile_yCenter = tile_yMin + (tile_yMax - tile_yMin) / 2;
    // determine whether the tile needs to be drawn
    short polygon_treeLevel = switch_level + 8;
    if (z <= polygon_treeLevel) {
        rtree_node = locRtreeNode(z, tile_xCenter, tile_yCenter, polygonRNode, switch_level);
    }
    else{
        rtree_node = locRtreeNode(polygon_treeLevel, tile_xCenter, tile_yCenter, polygonRNode, switch_level);
    }
    // non-blank tile
    if (rtree_node != nullptr){
        if (rtree_node->topo_type == 'i'){
            // get the features by spatial range retrieve
            vector<BOX_PAIR> polys;
            Box box(POINT(tile_xMin, tile_yMin), POINT(tile_xMax, tile_yMax));
            rtree_node->Rtree.query(bgi::intersects(box), back_inserter(polys));
            string data_string = "";
            for (int i = 0; i < polys.size(); i++){
                XY_VEC ring_cord = polys[i].second;
                data_string += "\"POLYGON((";
                for (int j = 0; j < ring_cord.size()-2; j+=2){
                    data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                }
                data_string += (to_string(ring_cord[0]) + " " + to_string(ring_cord[1]) + "))\"\n");
            }
            // plot features
            parameters params;
            params["type"] = "csv";
            params["inline"] = "wkt\n" + data_string;
            datasource_ptr polygon_ds = datasource_cache::instance().create(params);
            img_str = data_tileRender_a(tile_xMin, tile_yMin, tile_xMax, tile_yMax, R, G, B, AD, width, polygon_ds);
        }
        // fill tile
        else{
            Map map(256, 256);
            const color back_color(R, G, B, AD*255);
            map.set_background(back_color);
            image_rgba8 im(256, 256);
            agg_renderer<image_rgba8> ren(map, im);
            ren.apply();
            ostringstream im_os;
            save_to_stream(im, im_os, "png");
            img_str = im_os.str();
        }
        
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



double DiPGPlot(int z, int x, int y, float width, polygonNode *RNode){
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

    polygonNode *tile_node = locTileNode(z, x, y, RNode);
    if (tile_node->topo_type == 'i'){
        dis_tileRender(z, x, y, 255, 0, 0, 1, width, tile_node, RNode, tile_ptr);
    }
    else{
        # pragma omp parallel for num_threads(4)
        for (int i = 0; i < 256; i++){
            for (int j = 0; j < 256; j++){
                tile_ptr[i][4 * j] = 255;
                tile_ptr[i][4 * j + 1] = 0;
                tile_ptr[i][4 * j + 2] = 0;
                tile_ptr[i][4 * j + 3] = 255;
            }
        }
    }
    

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


double DaPGPlot(int z, int x, int y, float width, polygonNode *RNode, short transi_level){
    polygonNode *tile_node = locTileNode(z, x, y, RNode);

    if (tile_node->topo_type == 'i'){
        string data_string = "";
        double tile_Rz = L / pow(2, z-1);
        double tile_xMin = x * tile_Rz - L;
        double tile_xMax = (x + 1) * tile_Rz - L;
        double tile_yMin = L - (y + 1) * tile_Rz;
        double tile_yMax = L - y * tile_Rz;

        vector<polygonRtree> Rtree_list;
        tile_node->locRtree(Rtree_list, transi_level);

        Box box(POINT(tile_xMin, tile_yMin), POINT(tile_xMax, tile_yMax));

        for (auto rtree: Rtree_list){
            vector<BOX_PAIR> polys;
            rtree.query(bgi::within(box), back_inserter(polys));

            for (int i = 0; i < polys.size(); i++){
                XY_VEC ring_cord = polys[i].second;
                data_string += "\"POLYGON((";
                for (int j = 0; j < ring_cord.size()-2; j+=2){
                    data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                }
                data_string += (to_string(ring_cord[0]) + " " + to_string(ring_cord[1]) + "))\"\n");
            }
        }

        parameters params;
        params["type"] = "csv";
        params["inline"] = "wkt\n" + data_string;
        datasource_ptr polygon_ds = datasource_cache::instance().create(params);

        auto s_time = chrono::high_resolution_clock::now();

        // plot features
        string img_str = data_tileRender_a(tile_xMin, tile_yMin, tile_xMax, tile_yMax, 255, 0, 0, 1, width, polygon_ds);
    
        auto e_time = chrono::high_resolution_clock::now();

        chrono::duration<double> time_use = e_time - s_time;

        return time_use.count();
    }
    else{
        auto s_time = chrono::high_resolution_clock::now();

        Map map(256, 256);
        const color back_color(255, 0, 0, 255);
        map.set_background(back_color);
        image_rgba8 im(256, 256);
        agg_renderer<image_rgba8> ren(map, im);
        ren.apply();
        ostringstream im_os;
        save_to_stream(im, im_os, "png");
        string img_str = im_os.str();

        auto e_time = chrono::high_resolution_clock::now();

        chrono::duration<double> time_use = e_time - s_time;

        return time_use.count();
    }
}