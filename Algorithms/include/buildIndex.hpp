#include "hicore/core.hpp"
#include "hicore/core_creator.hpp"
#include "hicore/core_iterator.hpp"
#include "hicore/shp.hpp"
#include "indexNode.hpp"

# ifndef BUILDINDEX_HPP_
# define BUILDINDEX_HPP_

using namespace HiGIS::IO;
using namespace HiGIS::Core;

void coordTran(double &x, double &y);
short calTransiLevel(vector<string> shp_files, string data_espg);
pointNode* pointIndex(vector<string> shp_files, string data_espg, short transi_level);
linestringNode* linestringIndex(vector<string> shp_files, string data_espg, short transi_level);
polygonNode* polygonIndex(vector<string> shp_files, string data_espg, short transi_level);

# endif