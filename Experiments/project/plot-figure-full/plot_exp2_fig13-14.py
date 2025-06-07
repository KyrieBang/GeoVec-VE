import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import rgb2hex
import plotly.graph_objects as go
import plotly.io as pio
import seaborn as sns



def run_command(command):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        return result.stdout
    
    except subprocess.CalledProcessError as e:
        print(f"[command failed !]: {e}")
        print("[error info:]", e.stderr)


def extractInfo(out_info):
    info_list = out_info.strip().split('\n')

    time_list = [1000 * float(item) for item in info_list[:-1]]

    line = (info_list[-1]).split(" ")
    total_time = float(line[0])
    drawing_speed = float(line[1])

    return time_list, total_time, drawing_speed


def generateResult(data_path, shp_type):
    command = "mpirun -np 4 ../project/build/exp2_fig13-14" + " " + data_path + " " + shp_type
    out_info = run_command(command)
    time_list, total_time, drawing_speed = extractInfo(out_info)
   
    return time_list, total_time, drawing_speed


def plotFig13(fig_path, total_time_list, speed_list):
    datasets_name = ["P1", "P2", "P3", "P4", "L1", "L2", "L3", "L4", "A1", "A2", "A3"]
    
    fig, ax1 = plt.subplots(figsize=(12, 8))

    # plot bars
    ax1.grid(axis='y', linestyle='-', linewidth=2, color=(217/255, 217/255, 217/255), alpha=0.3)
    ax1.bar(datasets_name, total_time_list, color=(168/255, 196/255, 230/255), width=0.4, edgecolor=(89/255, 151/255, 208/255), linewidth=3, label='Time consumed to render 5000 tiles')
    ax1.tick_params(axis='x', labelcolor=(128/255, 128/255, 128/255), labelsize=17, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)
    ax1.set_ylabel('TIME (S)', color=(128/255, 128/255, 128/255), fontsize=18, fontweight='bold')
    ax1.tick_params(axis='y', labelcolor=(128/255, 128/255, 128/255), labelsize=17, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)
    
    for label in ax1.get_xticklabels():
        label.set_fontweight('bold')
    
    for label in ax1.get_yticklabels():
        label.set_fontweight('bold')

    ax1.legend(loc='upper center', bbox_to_anchor=(0.27, -0.07), frameon=False, ncol=2, prop={'size': 18, 'weight': 'bold'})

    # plot lines
    ax2 = ax1.twinx()
    ax2.plot(datasets_name, speed_list, color=(237/255, 125/255, 49/255), marker='^', markersize=10, linestyle="-", linewidth=4, label='Tile rendering speed')
    ax2.set_ylabel('SPEED (TILES/S)', color=(128/255, 128/255, 128/255), fontsize=18, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor=(128/255, 128/255, 128/255), labelsize=17, colors=rgb2hex((191/255, 191/255, 191/255)), width=3)

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

    ax2.legend(loc='upper center', bbox_to_anchor=(0.8, -0.07), frameon=False, ncol=2, prop={'size': 18, 'weight': 'bold'})

    
    plt.tight_layout()
    
    plt.savefig(fig_path, dpi=500)


def plotFig14(fig_path, all_time_df):
    plt.style.use('ggplot')
    fig, axs = plt.subplots()  

    Q3 = all_time_df.quantile(0.75)
    Q1 = all_time_df.quantile(0.25)
    IQR = Q3 - Q1
    Q_min = Q1 - 1.5 * IQR
    Q_max = Q3 + 1.5 * IQR
    cleaned_df = all_time_df[(all_time_df > Q_min) & (all_time_df < Q_max)]

    plt.ylabel('TIME (MS)', size=15)             
    plt.tick_params(labelsize=12, color="gray")

    fig = sns.violinplot(data=cleaned_df, common_norm=False, density_norm="width")

    fig.spines['left'].set_color("grey")
    fig.spines['right'].set_color("grey")
    fig.spines['top'].set_color("grey")
    fig.spines['bottom'].set_color("grey")

    plt.tight_layout()
    plt.savefig(fig_path, dpi = 1000)




if __name__ == "__main__":
    total_time_list = 11 * [0]
    speed_list = 11 * [0]
    
    # get the results of P1 dataset
    P1_path = "../../Datasets/point/P1"
    P1_type = "point"
    P1_time_list, P1_total_time, P1_speed = generateResult(P1_path, P1_type)
    total_time_list[0] = P1_total_time
    speed_list[0] = P1_speed

    # get the results of P2 dataset
    P2_path = "../../Datasets/point/P2"
    P2_type = "point"
    P2_time_list, P2_total_time, P2_speed = generateResult(P2_path, P2_type)
    total_time_list[1] = P2_total_time
    speed_list[1] = P2_speed

    # get the results of P3 dataset
    P3_path = "../../Datasets/point/P3"
    P3_type = "point"
    P3_time_list, P3_total_time, P3_speed = generateResult(P3_path, P3_type)
    total_time_list[2] = P3_total_time
    speed_list[2] = P3_speed

    # get the results of P4 dataset
    P4_path = "../../Datasets/point/P4"
    P4_type = "point"
    P4_time_list, P4_total_time, P4_speed = generateResult(P4_path, P4_type)
    total_time_list[3] = P4_total_time
    speed_list[3] = P4_speed

    # get the results of L1 dataset
    L1_path = "../../Datasets/linestring/L1"
    L1_type = "linestring"
    L1_time_list, L1_total_time, L1_speed = generateResult(L1_path, L1_type)
    total_time_list[4] = L1_total_time
    speed_list[4] = L1_speed

    # get the results of L2 dataset
    L2_path = "../../Datasets/linestring/L2"
    L2_type = "linestring"
    L2_time_list, L2_total_time, L2_speed = generateResult(L2_path, L2_type)
    total_time_list[5] = L2_total_time
    speed_list[5] = L2_speed

    # get the results of L3 dataset
    L3_path = "../../Datasets/linestring/L3"
    L3_type = "linestring"
    L3_time_list, L3_total_time, L3_speed = generateResult(L3_path, L3_type)
    total_time_list[6] = L3_total_time
    speed_list[6] = L3_speed

    # get the results of L4 dataset
    L4_path = "../../Datasets/linestring/L4"
    L4_type = "linestring"
    L4_time_list, L4_total_time, L4_speed = generateResult(L4_path, L4_type)
    total_time_list[7] = L4_total_time
    speed_list[7] = L4_speed

    # get the results of A1 dataset
    A1_path = "../../Datasets/polygon/A1"
    A1_type = "polygon"
    A1_time_list, A1_total_time, A1_speed = generateResult(A1_path, A1_type)
    total_time_list[8] = A1_total_time
    speed_list[8] = A1_speed

    # get the results of A2 dataset
    A2_path = "../../Datasets/polygon/A2"
    A2_type = "polygon"
    A2_time_list, A2_total_time, A2_speed = generateResult(A2_path, A2_type)
    total_time_list[9] = A2_total_time
    speed_list[9] = A2_speed

    # get the results of A3 dataset
    A3_path = "../../Datasets/polygon/A3"
    A3_type = "polygon"
    A3_time_list, A3_total_time, A3_speed = generateResult(A3_path, A3_type)
    total_time_list[10] = A3_total_time
    speed_list[10] = A3_speed



    # plot the results
    plotFig13("../outputs/full/fig13/fig13.png", total_time_list, speed_list)

    all_time_df = pd.DataFrame({'P1':P1_time_list, 'P2':P2_time_list, 'P3':P3_time_list, 'P4':P4_time_list, 'L1':L1_time_list, 'L2':L2_time_list, 'L3':L3_time_list, 'L4':L4_time_list, 'A1':A1_time_list, 'A2':A2_time_list, 'A3':A3_time_list})
    plotFig14("../outputs/full/fig14/fig14.png", all_time_df)