import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex



def run_command(command):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        return result.stdout
    
    except subprocess.CalledProcessError as e:
        print(f"[command failed !]: {e}")
        print("[error info:]", e.stderr)



def generateResult(data_path, shp_type):
    command = "mpirun -np 4 ../project/build/exp3_fig15" + " " + data_path + " " + shp_type
    out_info = run_command(command)
    
    info_list = out_info.strip().split('\n')
    total_time = float(info_list[0]) / 60.0
   
    return total_time


def generateGeoSparkViz(data_path, shp_type):
    command = "spark-submit" 
    command += (f" --conf spark.executor.cores=32")
    command += (" --conf spark.driver.memory=50g")
    command += (" --conf spark.executor.memory=50g")
    command += (" --conf spark.driver.maxResultSize=10g")
    command += " ../project/build/geosparkviz.jar"
    command += (" " + data_path + " " + shp_type + " fig15")
    out_info = run_command(command)

    info_list = out_info.strip().split('\n')
    total_time = float(info_list[0]) / 60.0

    return total_time


def plotFig15(fig_path, PATD_time_list, GSViz_time_list):
    datasets_name = ["P1", "P2", "P3", "P4", "L1", "L2", "L3", "L4", "A1", "A2", "A3"]
    
    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(len(datasets_name))
    width = 0.3

    # plot bars
    ax.grid(axis='y', linestyle='-', linewidth=2, color=(217/255, 217/255, 217/255), alpha=0.3)
    ax.bar(x-width/2, PATD_time_list, width, color=(168/255, 196/255, 230/255), edgecolor=(89/255, 151/255, 208/255), linewidth=3, label='PATD')
    ax.bar(x+width/2, GSViz_time_list, width, color=(246/255, 172/255, 141/255), edgecolor=(232/255, 122/255, 48/255), linewidth=3, label='GeoSparkViz')
    ax.tick_params(axis='x', labelcolor=(128/255, 128/255, 128/255), labelsize=17, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)
    ax.set_ylabel('TIME (MIN)', color=(128/255, 128/255, 128/255), fontsize=18, fontweight='bold')
    ax.tick_params(axis='y', labelcolor=(128/255, 128/255, 128/255), labelsize=17, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)
    
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')

    ax.set_xticks(x)
    ax.set_xticklabels(datasets_name)
    legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07), frameon=False, ncol=3, prop={'size': 18, 'weight': 'bold'})
    for text in legend.get_texts():
        text.set_color((128/255, 128/255, 128/255))

    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax.spines['top'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax.spines['left'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)

    
    plt.tight_layout()
    
    plt.savefig(fig_path, dpi=500)







if __name__ == "__main__":
    PATD_time_list = 11 * [0]
    GeoSparkViz_time_list = 11 * [0]
    
    # get the results of P1 dataset
    P1_path = "../../Datasets/point/P1"
    P1_type = "point"
    P1_total_time = generateResult(P1_path, P1_type)
    PATD_time_list[0] = P1_total_time
    P1_total_time = generateGeoSparkViz(P1_path, P1_type)
    GeoSparkViz_time_list[0] = P1_total_time

    # get the results of P2 dataset
    P2_path = "../../Datasets/point/P2"
    P2_type = "point"
    P2_total_time = generateResult(P2_path, P2_type)
    PATD_time_list[1] = P2_total_time
    P2_total_time = generateGeoSparkViz(P2_path, P2_type)
    GeoSparkViz_time_list[1] = P2_total_time

    # get the results of P3 dataset
    P3_path = "../../Datasets/point/P3"
    P3_type = "point"
    P3_total_time = generateResult(P3_path, P3_type)
    PATD_time_list[2] = P3_total_time
    P3_total_time = generateGeoSparkViz(P3_path, P3_type)
    GeoSparkViz_time_list[2] = P3_total_time

    # get the results of P4 dataset
    P4_path = "../../Datasets/point/P4"
    P4_type = "point"
    P4_total_time = generateResult(P4_path, P4_type)
    PATD_time_list[3] = P4_total_time
    P4_total_time = generateGeoSparkViz(P4_path, P4_type)
    GeoSparkViz_time_list[3] = P4_total_time

    # get the results of L1 dataset
    L1_path = "../../Datasets/linestring/L1"
    L1_type = "linestring"
    L1_total_time = generateResult(L1_path, L1_type)
    PATD_time_list[4] = L1_total_time
    L1_total_time = generateGeoSparkViz(L1_path, L1_type)
    GeoSparkViz_time_list[4] = L1_total_time

    # get the results of L2 dataset
    L2_path = "../../Datasets/linestring/L2"
    L2_type = "linestring"
    L2_total_time = generateResult(L2_path, L2_type)
    PATD_time_list[5] = L2_total_time
    L2_total_time = generateGeoSparkViz(L2_path, L2_type)
    GeoSparkViz_time_list[5] = L2_total_time

    # get the results of L3 dataset
    L3_path = "../../Datasets/linestring/L3"
    L3_type = "linestring"
    L3_total_time = generateResult(L3_path, L3_type)
    PATD_time_list[6] = L3_total_time
    L3_total_time = generateGeoSparkViz(L3_path, L3_type)
    GeoSparkViz_time_list[6] = L3_total_time

    # get the results of L4 dataset
    L4_path = "../../Datasets/linestring/L4"
    L4_type = "linestring"
    L4_total_time = generateResult(L4_path, L4_type)
    PATD_time_list[7] = L4_total_time
    L4_total_time = generateGeoSparkViz(L4_path, L4_type)
    GeoSparkViz_time_list[7] = L4_total_time

    # get the results of A1 dataset
    A1_path = "../../Datasets/polygon/A1"
    A1_type = "polygon"
    A1_total_time = generateResult(A1_path, A1_type)
    PATD_time_list[8] = A1_total_time
    A1_total_time = generateGeoSparkViz(A1_path, A1_type)
    GeoSparkViz_time_list[8] = A1_total_time

    # get the results of A2 dataset
    A2_path = "../../Datasets/polygon/A2"
    A2_type = "polygon"
    A2_total_time = generateResult(A2_path, A2_type)
    PATD_time_list[9] = A2_total_time
    A2_total_time = generateGeoSparkViz(A2_path, A2_type)
    GeoSparkViz_time_list[9] = A2_total_time

    # get the results of A3 dataset
    A3_path = "../../Datasets/polygon/A3"
    A3_type = "polygon"
    A3_total_time = generateResult(A3_path, A3_type)
    PATD_time_list[10] = A3_total_time
    A3_total_time = generateGeoSparkViz(A3_path, A3_type)
    GeoSparkViz_time_list[10] = A3_total_time

    
    # plot the results
    plotFig15("../outputs/full/fig15/fig15.png", PATD_time_list, GeoSparkViz_time_list)