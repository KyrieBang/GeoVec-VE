import subprocess
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex


def run_command(command):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        return result.stdout
    
    except subprocess.CalledProcessError as e:
        print(f"[command failed !]: {e}")
        print("[error info:]", e.stderr)


def extractInfo(out_info):
    info_list = out_info.strip().split('\n')

    speed_list = []

    for speed in info_list:
       speed_list.append(float(speed))

    return speed_list


def extractInfo_GV(out_info):
    info_list = out_info.strip().split('\n')

    speed_list = []

    for speed in info_list:
       speed_list.append(float(speed))

    return speed_list


def plotSubFig(tile_label, fig_path, PATD_speed_list, GSViz_speed_list):

    x_list = ["level 1", "level 3", "level 5", "level 7", "level 9"]
    
    fig, ax = plt.subplots()

    # plot lines
    ax.tick_params(axis='x', labelcolor=(128/255, 128/255, 128/255), labelsize=13, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)
    ax.set_ylabel('Drawing Speed (tiles/s)', color=(128/255, 128/255, 128/255), fontsize=13, fontweight='bold')
    ax.tick_params(axis='y', labelcolor=(128/255, 128/255, 128/255), labelsize=13, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)

    ax.grid(axis='y', linestyle='-', linewidth=2, color=(217/255, 217/255, 217/255), alpha=0.3)
    ax.plot(x_list, PATD_speed_list, color=(91/255, 155/255, 213/255), marker='o', markersize=10, linestyle="-", linewidth=4, label='PATD')
    ax.plot(x_list, GSViz_speed_list, color=(237/255, 125/255, 49/255), marker='^', markersize=10, linestyle="-", linewidth=4, label='GeoSparkViz')
    

    
   
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')


    ax.spines['bottom'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax.spines['left'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    
    fig.suptitle(tile_label, fontsize=15, color=(0, 0, 0), fontweight='bold', y=0.05)
    fig.subplots_adjust(top=0.95, bottom=0.2, left=0.13, right=0.97)
     
    plt.savefig(fig_path, dpi=500)


def generateResult(data_path, shp_type, fig_path, tile_label):
    command = "mpirun -np 2 ../project/build/exp3_fig16" + " " + data_path + " " + shp_type
    out_info = run_command(command)
    PATD_speed_list = extractInfo(out_info)

    command2 = "spark-submit" 
    command2 += (" --conf spark.executor.cores=16")
    command2 += (" --conf spark.driver.memory=50g")
    command2 += (" --conf spark.executor.memory=50g")
    command2 += (" --conf spark.driver.maxResultSize=10g")
    command2 += " ../project/build/geosparkviz.jar"
    command2 += (" " + data_path + " " + shp_type + " fig16")
    out_info2 = run_command(command2)
    GSViz_speed_list = extractInfo_GV(out_info2)
    
    plotSubFig(tile_label, fig_path, PATD_speed_list, GSViz_speed_list)




if __name__ == "__main__":
    # get the results of P4 dataset
    P4_path = "../../Datasets/point/P4"
    P4_type = "point"
    P4_fig_path = "../outputs/simple/fig16/P4.png"
    generateResult(P4_path, P4_type, P4_fig_path, "(a) $P_4$")

    # get the results of L2 dataset
    L2_path = "../../Datasets/linestring/L2"
    L2_type = "linestring"
    L2_fig_path = "../outputs/simple/fig16/L2.png"
    generateResult(L2_path, L2_type, L2_fig_path, "(b) $L_2$")
    
    # get the results of A2 dataset
    A2_path = "../../Datasets/polygon/A2"
    A2_type = "polygon"
    A2_fig_path = "../outputs/simple/fig16/A2.png"
    generateResult(A2_path, A2_type, A2_fig_path, "(e) $A_2$")
