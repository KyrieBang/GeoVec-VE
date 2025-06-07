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

    DiPG_time_list = []
    DaPG_time_list = []
    DiPG_speed_list = []
    DaPG_speed_list = []

    for s_info in info_list:
       line = s_info.split(" ")
       DiPG_time_list.append(float(line[0]))
       DaPG_time_list.append(float(line[1]))
       DiPG_speed_list.append(float(line[2]))
       DaPG_speed_list.append(float(line[3]))

    return DiPG_time_list, DaPG_time_list, DiPG_speed_list, DaPG_speed_list


def plotSubFig(tile_label, min_level, max_level, fig_path, DiPG_time_list, DaPG_time_list, DiPG_speed_list, DaPG_speed_list):
    zoom_levels = [str(i) for i in range(min_level, max_level + 1)]
    
    fig, ax1 = plt.subplots()

    # plot bars
    x = np.arange(len(zoom_levels))
    width = 0.3

    ax1.grid(axis='y', linestyle='-', linewidth=2, color=(217/255, 217/255, 217/255), alpha=0.3)
    ax1.bar(x - width/2, DiPG_time_list, width, color=(168/255, 196/255, 230/255), edgecolor=(89/255, 151/255, 208/255), linewidth=3, label='Total drawing time of DiPG')
    ax1.bar(x + width/2, DaPG_time_list, width, color=(246/255, 172/255, 141/255), edgecolor=(232/255, 122/255, 48/255), linewidth=3, label='Total drawing time of DaPG')
    ax1.set_xlabel('ZOOM LEVEL', color=(128/255, 128/255, 128/255), fontsize=14, fontweight='bold')
    ax1.tick_params(axis='x', labelcolor=(128/255, 128/255, 128/255), labelsize=13, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)
    ax1.set_ylabel('TIME (S)', color=(128/255, 128/255, 128/255), fontsize=14, fontweight='bold')
    ax1.tick_params(axis='y', labelcolor=(128/255, 128/255, 128/255), labelsize=13, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)
    
    for label in ax1.get_xticklabels():
        label.set_fontweight('bold')
    
    for label in ax1.get_yticklabels():
        label.set_fontweight('bold')

    # plot lines
    ax2 = ax1.twinx()
    ax2.plot(zoom_levels, DiPG_speed_list, color=(91/255, 155/255, 213/255), marker='o', markersize=10, linestyle=":", linewidth=4, label='Tile drawing speed of DiPG')
    ax2.plot(zoom_levels, DaPG_speed_list, color=(237/255, 125/255, 49/255), marker='^', markersize=10, linestyle=":", linewidth=4, label='Tile drawing speed of DiPG')
    ax2.set_ylabel('SPEED (TILES/S)', color=(128/255, 128/255, 128/255), fontsize=13, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor=(128/255, 128/255, 128/255), labelsize=13, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)

    for label in ax2.get_yticklabels():
        label.set_fontweight('bold')

    ax2.spines['bottom'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax2.spines['top'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax2.spines['left'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax2.spines['right'].set_color(rgb2hex((191/255, 191/255, 191/255)))
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['top'].set_linewidth(2)
    ax2.spines['left'].set_linewidth(2)
    ax2.spines['right'].set_linewidth(2)
    
    fig.suptitle(tile_label, fontsize=15, color=(0, 0, 0), fontweight='bold', y=0.05)
    fig.subplots_adjust(top=0.95, bottom=0.2, left=0.13, right=0.87)
     
    plt.savefig(fig_path, dpi=500)


def generateResult(data_path, shp_type, fig_path, tile_label, min_level, max_level):
    command = "../project/build/exp1_fig8" + " " + data_path + " " + shp_type + " " + str(min_level) + " " + str(max_level)
    out_info = run_command(command)
    DiPG_time_list, DaPG_time_list, DiPG_speed_list, DaPG_speed_list = extractInfo(out_info)
    plotSubFig(tile_label, min_level, max_level, fig_path, DiPG_time_list, DaPG_time_list, DiPG_speed_list, DaPG_speed_list)




if __name__ == "__main__":
    # get the subfig results of P1 dataset in fig.8
    P1_path = "../../Datasets/point/P1"
    P1_type = "point"
    P1_fig_path = "../outputs/full/fig8/P1.png"
    P1_min_level = 11
    P1_max_level = 16
    generateResult(P1_path, P1_type, P1_fig_path, "(a) $P_1$", P1_min_level, P1_max_level)

    # get the subfig results of P2 dataset in fig.8
    P2_path = "../../Datasets/point/P2"
    P2_type = "point"
    P2_fig_path = "../outputs/full/fig8/P2.png"
    P2_min_level = 10
    P2_max_level = 15
    generateResult(P2_path, P2_type, P2_fig_path, "(b) $P_2$", P2_min_level, P2_max_level)

    # get the subfig results of P3 dataset in fig.8
    P3_path = "../../Datasets/point/P3"
    P3_type = "point"
    P3_fig_path = "../outputs/full/fig8/P3.png"
    P3_min_level = 7
    P3_max_level = 12
    generateResult(P3_path, P3_type, P3_fig_path, "(c) $P_3$", P3_min_level, P3_max_level)

    # get the subfig results of P4 dataset in fig.8
    P4_path = "../../Datasets/point/P4"
    P4_type = "point"
    P4_fig_path = "../outputs/full/fig8/P4.png"
    P4_min_level = 6
    P4_max_level = 11
    generateResult(P4_path, P4_type, P4_fig_path, "(d) $P_4$", P4_min_level, P4_max_level)

    # get the subfig results of L1 dataset in fig.8
    L1_path = "../../Datasets/linestring/L1"
    L1_type = "linestring"
    L1_fig_path = "../outputs/full/fig8/L1.png"
    L1_min_level = 11
    L1_max_level = 16
    generateResult(L1_path, L1_type, L1_fig_path, "(e) $L_1$", L1_min_level, L1_max_level)

    # get the subfig results of L2 dataset in fig.8
    L2_path = "../../Datasets/linestring/L2"
    L2_type = "linestring"
    L2_fig_path = "../outputs/full/fig8/L2.png"
    L2_min_level = 7
    L2_max_level = 12
    generateResult(L2_path, L2_type, L2_fig_path, "(f) $L_2$", L2_min_level, L2_max_level)

    # get the subfig results of L3 dataset in fig.8
    L3_path = "../../Datasets/linestring/L3"
    L3_type = "linestring"
    L3_fig_path = "../outputs/full/fig8/L3.png"
    L3_min_level = 7
    L3_max_level = 12
    generateResult(L3_path, L3_type, L3_fig_path, "(g) $L_3$", L3_min_level, L3_max_level)

    # get the subfig results of A1 dataset in fig.8
    A1_path = "../../Datasets/polygon/A1"
    A1_type = "polygon"
    A1_fig_path = "../outputs/full/fig8/A1.png"
    A1_min_level = 10
    A1_max_level = 15
    generateResult(A1_path, A1_type, A1_fig_path, "(h) $A_1$", A1_min_level, A1_max_level)

    # get the subfig results of A2 dataset in fig.8
    A2_path = "../../Datasets/polygon/A2"
    A2_type = "polygon"
    A2_fig_path = "../outputs/full/fig8/A2.png"
    A2_min_level = 8
    A2_max_level = 13
    generateResult(A2_path, A2_type, A2_fig_path, "(i) $A_2$", A2_min_level, A2_max_level)