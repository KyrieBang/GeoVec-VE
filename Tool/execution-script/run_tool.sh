# !/bin/bash
point_shpPath="../../Datasets/point/POI"
linestring_shpPath="../../Datasets/linestring/Railway"
polygon_shpPath="../../Datasets/polygon/Building"

nohup ../visual-engine/visualEngine $point_shpPath $linestring_shpPath $polygon_shpPath > ./tool_message.log 2>&1 &
python3 ../browser-interface/run_interface.py