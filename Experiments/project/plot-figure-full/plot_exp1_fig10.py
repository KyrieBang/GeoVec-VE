import subprocess
import matplotlib.pyplot as plt
import numpy as np
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

    DiPG_speed_list = []
    DaPG_speed_list = []

    for s_info in info_list:
       line = s_info.split(" ")
       DiPG_speed_list.append(float(line[0]))
       DaPG_speed_list.append(float(line[1]))

    return DiPG_speed_list, DaPG_speed_list


def plotSubFig(tile_label, min_len, max_len, fig_path, DiPG_speed_list, DaPG_speed_list):
    num_list = []
    num_list.append("0-" + str(int(min_len/100)))
    for i in range(int(min_len/100), int(max_len/100)):
        num_list.append(str(i) + "-" + str(i+1))
    num_list.append(">" + str(int(max_len/100)))
    
    fig, ax = plt.subplots()

    # plot lines
    ax.set_xlabel("f_len ($10^2$KM)", color=(128/255, 128/255, 128/255), fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', labelcolor=(128/255, 128/255, 128/255), labelsize=13, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)
    ax.set_ylabel('Drawing Speed (tiles/s)', color=(128/255, 128/255, 128/255), fontsize=13, fontweight='bold')
    ax.tick_params(axis='y', labelcolor=(128/255, 128/255, 128/255), labelsize=13, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)

    ax.grid(axis='y', linestyle='-', linewidth=2, color=(217/255, 217/255, 217/255), alpha=0.3)
    ax.plot(num_list, DiPG_speed_list, color=(91/255, 155/255, 213/255), marker='o', markersize=10, linestyle="-", linewidth=4, label='tile drawing speed of DiPG')
    ax.plot(num_list, DaPG_speed_list, color=(237/255, 125/255, 49/255), marker='^', markersize=10, linestyle="-", linewidth=4, label='Tile drawing speed of DiPG')
    

    
   
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


def generateResult(data_path, shp_type, fig_path, tile_label, min_len, max_len, min_level, max_level):
    command = "mpirun -np 4 ../project/build/exp1_fig10" + " " + data_path + " " + shp_type + " " + str(min_len) + " " + str(max_len) + " " + str(min_level) + " " + str(max_level)
    out_info = run_command(command)
    DiPG_speed_list, DaPG_speed_list = extractInfo(out_info)
    
    plotSubFig(tile_label, min_len, max_len, fig_path, DiPG_speed_list, DaPG_speed_list)




if __name__ == "__main__":
    # get the subfig results of L1 dataset in fig.10
    L1_path = "../../Datasets/linestring/L1"
    L1_type = "linestring"
    L1_fig_path = "../outputs/fig10/L1.png"
    L1_min_len = 200
    L1_max_len = 600
    L1_min_level = 11
    L1_max_level = 16
    generateResult(L1_path, L1_type, L1_fig_path, "(a) $L_1$", L1_min_len, L1_max_len, L1_min_level, L1_max_level)

    # get the subfig results of L2 dataset in fig.10
    L2_path = "../../Datasets/linestring/L2"
    L2_type = "linestring"
    L2_fig_path = "../outputs/fig10/L2.png"
    L2_min_len = 200
    L2_max_len = 600
    L2_min_level = 7
    L2_max_level = 12
    generateResult(L2_path, L2_type, L2_fig_path, "(b) $L_2$", L2_min_len, L2_max_len, L2_min_level, L2_max_level)

    # get the subfig results of L3 dataset in fig.10
    L3_path = "../../Datasets/linestring/L3"
    L3_type = "linestring"
    L3_fig_path = "../outputs/fig10/L3.png"
    L3_min_len = 200
    L3_max_len = 600
    L3_min_level = 7
    L3_max_level = 10
    generateResult(L3_path, L3_type, L3_fig_path, "(c) $L_3$", L3_min_len, L3_max_len, L3_min_level, L3_max_level)

    # get the subfig results of A1 dataset in fig.10
    A1_path = "../../Datasets/polygon/A1"
    A1_type = "polygon"
    A1_fig_path = "../outputs/fig10/A1.png"
    A1_min_len = 300
    A1_max_len = 700
    A1_min_level = 10
    A1_max_level = 15
    generateResult(A1_path, A1_type, A1_fig_path, "(d) $A_1$", A1_min_len, A1_max_len, A1_min_level, A1_max_level)

    # get the subfig results of A2 dataset in fig.10
    A2_path = "../../Datasets/polygon/A2"
    A2_type = "polygon"
    A2_fig_path = "../outputs/fig10/A2.png"
    A2_min_len = 300
    A2_max_len = 700
    A2_min_level = 8
    A2_max_level = 13
    generateResult(A2_path, A2_type, A2_fig_path, "(e) $A_2$", A2_min_len, A2_max_len, A2_min_level, A2_max_level)