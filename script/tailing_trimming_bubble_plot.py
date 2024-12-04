#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :tailing_trimming_bubble_plot.py
# @Time      :2023/12/20 20:01:22
# @Author    :Yuchen@rlab
# @Description: This script is used to generate bubble plot for tailing trimming result
# @version   :1

import yaml
import sys
import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerPathCollection, HandlerLine2D
from multiprocessing import Process

# 忽略RuntimeWarning类型的警告
warnings.filterwarnings("ignore", category=RuntimeWarning)

def plot_matrix(df_list, total_list, num_cols):
    num_rows = len(df_list) // num_cols

    # 创建一个新的图形，并设置子图布局
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 10))
    axes = np.reshape(axes, (num_rows, num_cols))
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'cyan']
    # 遍历每个DataFrame
    for i, df in enumerate(df_list):
        # 删除最后一列
        df = df.drop(df.columns[:3], axis=1).drop(df.columns[-1], axis=1)
        df = df.drop(df.index[:3])

        # 将DataFrame转换为矩阵
        matrix = df.values
        matrix_sum = np.sum(matrix)

        # 获取矩阵的行数和列数
        rows, cols = matrix.shape

        # 计算每个点的坐标
        x = np.arange(cols)
        y = np.arange(rows)
        X, Y = np.meshgrid(x, y)

        # 将矩阵展平为一维数组，并计算点的半径
        values = matrix.flatten()
        radii = (values / matrix_sum) * 1500

        # 倒转y轴坐标
        axes[i // num_cols, i % num_cols].invert_yaxis()
        axes[i // num_cols, i % num_cols].yaxis.set_label_position('right') 

        x_ticks = np.arange(8)   # 生成从0到10的数组
        x_labels = np.arange(7, -1, -1) # 生成从10到0的数组
        axes[i // num_cols, i % num_cols].set_xticks(x_ticks)
        axes[i // num_cols, i % num_cols].set_xticklabels(x_labels)

        y_ticks = np.arange(8)  # 生成从0到10的数组
        y_labels = np.arange(7, -1, -1) # 生成从0到10的数组
        axes[i // num_cols, i % num_cols].set_yticks(y_ticks)
        axes[i // num_cols, i % num_cols].set_yticklabels(y_labels)

        ax = axes[i // num_cols, i % num_cols]
        ax.yaxis.tick_right()
        ax.tick_params(axis='both', which='both', size=0)
        ax.set_aspect('equal')
        ax.margins(x=0.1, y=0.1)

        axes[i // num_cols, i % num_cols].grid(True, zorder=0)
        axes[i // num_cols, i % num_cols].scatter(X.flatten(), Y.flatten(), s=radii, zorder=1000, color=colors[i % len(colors)])
        axes[i // num_cols, i % num_cols].set_title(f'{df_list[i].name} ({total_list[i]})')
        axes[i // num_cols, i % num_cols].set_xlabel('Length of trimming',labelpad=10,fontsize=12)
        axes[i // num_cols, i % num_cols].set_ylabel('Length of tailing',labelpad=10,rotation=-90,fontsize=12)

    # 调整子图布局
    # fig.text(0.5,0.2, 'Length of trimming', ha='center',fontsize=12)
    # fig.text(1.01,0.5, 'Length of tailing', va='center', rotation=-90, fontsize=12)
    fig.suptitle(mir_name,fontsize=12)
    fig.tight_layout()
    plt.savefig(os.path.join(plot_output, mir_name + '_matrix_plot.pdf'), format='pdf')
    # plt.show()
    plt.close()

def plot_matrix_merge(df_list):
    # 创建一个新的图形
    plt.figure(figsize=(10, 6))
    for i, df in enumerate(df_list):
        # 删除最后一列
        df = df.drop(df.columns[:3], axis=1).drop(df.columns[-1], axis=1)
        df = df.drop(df.index[:3])

        # 将DataFrame转换为矩阵
        matrix = df.values
        matrix_sum = np.sum(matrix)

        # 获取矩阵的行数和列数
        rows, cols = matrix.shape
        
        # 计算每个点的坐标
        x = np.arange(cols)
        y = np.arange(rows)
        X, Y = np.meshgrid(x, y)

        # 将矩阵展平为一维数组，并计算点的半径
        values = matrix.flatten()
        radii = (values / matrix_sum) * 1500

        
        plt.gca().yaxis.set_label_position('right')
        
        x_ticks = np.arange(8)   # 生成从0到7的数组
        x_labels = np.arange(7, -1, -1) # 生成从10到0的数组
        plt.xticks(x_ticks, x_labels)

        y_ticks = np.arange(8)  # 生成从0到10的数组
        y_labels = np.arange(7, -1, -1) # 生成从0到10的数组
        plt.yticks(y_ticks, y_labels)

        ax = plt.gca()
        ax.yaxis.tick_right()
        ax.set_aspect('equal')
        ax.margins(x=0.1, y=0.1)

        plt.grid(True, zorder=0)
        plt.scatter(X.flatten(), Y.flatten(), s=radii, zorder=1000, label=f'{df_list[i].name}\n({total_list[i]})',alpha=0.5) 
        plt.xlabel('Length of trimming',labelpad=10,fontsize=12)
        plt.ylabel('Length of tailing',labelpad=10,rotation=-90,fontsize=12)
        plt.title(mir_name,fontsize=12)
        plt.legend(handler_map={PathCollection: HandlerPathCollection(update_func=update_scatter),plt.Line2D: HandlerLine2D(update_func=updateline)},loc='lower center', bbox_to_anchor=(0.5, -0.3),ncol=len(df_list), prop={'size': 9})
    
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_output, mir_name + '_matrix_plot_merge.pdf'), format='pdf')
    # plt.show()
    plt.close()

def update_scatter(handle, orig):
    handle.update_from(orig)
    handle.set_sizes([64])

def updateline(handle, orig):
    handle.update_from(orig)
    handle.set_markersize(8)

if __name__ == '__main__':
# 从配置文件中加载配置
    with open(sys.argv[1], 'r') as config_file:
        config = yaml.safe_load(config_file)

    mir_list = pd.read_csv(config['meta_file'], sep='\t', header=None)
    num_cols = config['num_cols']
    plot_output = config['plot_output']

    progress1 = []
    progress2 = []
    for mir in mir_list[0]:
        mir_name = mir.replace('*', '_star')
        df_list = []
        total_list = []
        for df_config in config['df']:
            df_path = df_config['path']
            df_name = df_config['name']
            matrix_file = os.path.join(df_path, mir_name + '.txt')
            df = pd.read_csv(matrix_file, sep='\t', skiprows=1, header=None)
            df.name = df_name
            df_list.append(df)

            summary = pd.read_csv(matrix_file, sep='\t', nrows=1, header=None)
            total_count = summary[1][0]
            total_list.append(total_count)
        plot_matrix(df_list, total_list,num_cols=num_cols)
        plot_matrix_merge(df_list)

        p1 = Process(target=plot_matrix, args=(df_list, total_list,num_cols))
        p1.start()
        progress1.append(p1)
        p2 = Process(target=plot_matrix_merge, args=(df_list,))
        p2.start()
        progress2.append(p2)
    for p1 in progress1:
        p1.join()
    for p2 in progress2:
        p2.join()

    print("All done!")