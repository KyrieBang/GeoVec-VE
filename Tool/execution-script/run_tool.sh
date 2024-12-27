# !/bin/bash
point_shpPath="../../Datasets/point/China_point"
linestring_shpPath="../../Datasets/linestring/China_railway"
polygon_shpPath="../../Datasets/polygon/China_building"

nohup ../visual-engine/visualEngine $point_shpPath $linestring_shpPath $polygon_shpPath > ./tool_message.log 2>&1 &
python3 ../browser-interface/run_interface.py