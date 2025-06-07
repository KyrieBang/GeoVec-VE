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


def plotSubFig(tile_label, min_count, max_count, fig_path, DiPG_speed_list, DaPG_speed_list):
    num_list = []
    num_list.append("0-" + str(int(min_count/100)))
    for i in range(int(min_count/100), int(max_count/100)):
        num_list.append(str(i) + "-" + str(i+1))
    num_list.append(">" + str(int(max_count/100)))
    
    fig, ax = plt.subplots()

    # plot lines
    ax.set_xlabel("f_count ($10^2$)", color=(128/255, 128/255, 128/255), fontsize=14, fontweight='bold')
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


def generateResult(data_path, shp_type, fig_path, tile_label, min_count, max_count, min_level, max_level):
    command = "mpirun -np 4 ../project/build/exp1_fig9" + " " + data_path + " " + shp_type + " " + str(min_count) + " " + str(max_count) + " " + str(min_level) + " " + str(max_level)
    out_info = run_command(command)
    DiPG_speed_list, DaPG_speed_list = extractInfo(out_info)
    
    plotSubFig(tile_label, min_count, max_count, fig_path, DiPG_speed_list, DaPG_speed_list)




if __name__ == "__main__":
    # get the subfig results of P1 dataset in fig.9
    P1_path = "../../Datasets/point/P1"
    P1_type = "point"
    P1_fig_path = "../outputs/simple/fig9/P1.png"
    P1_min_count = 600
    P1_max_count = 1000
    P1_min_level = 11
    P1_max_level = 16
    generateResult(P1_path, P1_type, P1_fig_path, "(a) $P_1$", P1_min_count, P1_max_count, P1_min_level, P1_max_level)

    # get the subfig results of P2 dataset in fig.9
    P2_path = "../../Datasets/point/P2"
    P2_type = "point"
    P2_fig_path = "../outputs/simple/fig9/P2.png"
    P2_min_count = 600
    P2_max_count = 1000
    P2_min_level = 10
    P2_max_level = 15
    generateResult(P2_path, P2_type, P2_fig_path, "(b) $P_2$", P2_min_count, P2_max_count, P2_min_level, P2_max_level)

    # get the subfig results of P3 dataset in fig.9
    P3_path = "../../Datasets/point/P3"
    P3_type = "point"
    P3_fig_path = "../outputs/simple/fig9/P3.png"
    P3_min_count = 600
    P3_max_count = 1000
    P3_min_level = 7
    P3_max_level = 12
    generateResult(P3_path, P3_type, P3_fig_path, "(c) $P_3$", P3_min_count, P3_max_count, P3_min_level, P3_max_level)

    # get the subfig results of P4 dataset in fig.9
    P4_path = "../../Datasets/point/P4"
    P4_type = "point"
    P4_fig_path = "../outputs/simple/fig9/P4.png"
    P4_min_count = 600
    P4_max_count = 1000
    P4_min_level = 6
    P4_max_level = 11
    generateResult(P4_path, P4_type, P4_fig_path, "(d) $P_4$", P4_min_count, P4_max_count, P4_min_level, P4_max_level)

    # get the subfig results of L1 dataset in fig.9
    L1_path = "../../Datasets/linestring/L1"
    L1_type = "linestring"
    L1_fig_path = "../outputs/simple/fig9/L1.png"
    L1_min_count = 600
    L1_max_count = 1000
    L1_min_level = 11
    L1_max_level = 16
    generateResult(L1_path, L1_type, L1_fig_path, "(e) $L_1$", L1_min_count, L1_max_count, L1_min_level, L1_max_level)

    # get the subfig results of L2 dataset in fig.9
    L2_path = "../../Datasets/linestring/L2"
    L2_type = "linestring"
    L2_fig_path = "../outputs/simple/fig9/L2.png"
    L2_min_count = 600
    L2_max_count = 1000
    L2_min_level = 7
    L2_max_level = 12
    generateResult(L2_path, L2_type, L2_fig_path, "(f) $L_2$", L2_min_count, L2_max_count, L2_min_level, L2_max_level)

    # get the subfig results of A1 dataset in fig.9
    A1_path = "../../Datasets/polygon/A1"
    A1_type = "polygon"
    A1_fig_path = "../outputs/simple/fig9/A1.png"
    A1_min_count = 1800
    A1_max_count = 2200
    A1_min_level = 10
    A1_max_level = 15
    generateResult(A1_path, A1_type, A1_fig_path, "(h) $A_1$", A1_min_count, A1_max_count, A1_min_level, A1_max_level)

    # get the subfig results of A2 dataset in fig.9
    A2_path = "../../Datasets/polygon/A2"
    A2_type = "polygon"
    A2_fig_path = "../outputs/simple/fig9/A2.png"
    A2_min_count = 1800
    A2_max_count = 2200
    A2_min_level = 8
    A2_max_level = 13
    generateResult(A2_path, A2_type, A2_fig_path, "(i) $A_2$", A2_min_count, A2_max_count, A2_min_level, A2_max_level)
    