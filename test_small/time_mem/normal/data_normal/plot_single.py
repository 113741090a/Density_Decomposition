import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
from matplotlib.table import Table
import os

# 读取数据文件的函数
def read_data(filename):
    sizes = []
    averages = []

    with open(filename, 'r') as file:
        for line in file:
            data = line.split()
            if len(data) != 4:
                continue  # 忽略格式错误的行
            average = float(data[0])/100000
            size = int(data[3])  # 第四列是 size

            averages.append(average)
            sizes.append(size)

    return sizes, averages

# 定义文件的颜色映射
file_colors = {
    'PR_lin': 'magenta',
    'PR': 'orange',
    'PR_exp': 'red',
    'FW_Elist': 'lightblue',
    'Elist++': 'green',
    'Greedy++': 'darkblue',
    'FISTA*': 'cyan'
}

def plot_data(datasets, output_dir):
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('tight')
    ax.axis('off')
    
    # 准备表格数据
    table_data = []
    
    # 表头 - 使用实际的文件大小作为列名
    if datasets:
        sizes = datasets[0][0]  # 使用第一个数据集的sizes作为列名
        headers = ['Algorithm'] + [f'Size |F| {size}' for size in sizes]
    else:
        headers = ['Algorithm']
    
    table_data.append(headers)
    
    # 添加每个数据集的数据
    for sizes, averages, name, color in datasets:
        row = [name] + [f'{avg:.4f}s' for avg in averages]  # 时间单位设为秒
        table_data.append(row)
    
    # 创建表格
    table = ax.table(cellText=table_data, 
                    loc='center', 
                    cellLoc='center')
    
    # 设置表格样式
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # 设置表头样式
    for i in range(len(headers)):
        table[(0, i)].set_facecolor('#DDDDDD')
        table[(0, i)].set_text_props(weight='bold')
    
    # 为不同算法设置行颜色
    # for i, (_, _, name, color) in enumerate(datasets, 1):
    #     for j in range(len(headers)):
    #         table[(i, j)].set_facecolor(color)
    
    # 设置标题
    plt.title('Algorithm Time Comparison (seconds)', fontsize=14, fontweight='bold', pad=20)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存表格
    output_path = os.path.join(output_dir, 'time_comparison_table.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    #plt.show()
    
    print(f"时间比较表格已保存到: {output_path}")
    
    # 同时创建DataFrame并保存为CSV
    # df_data = {}
    # for sizes, averages, name, color in datasets:
    #     df_data[name] = [f'{avg:.4f}s' for avg in averages]
    
    # if datasets:
    #     df = pd.DataFrame(df_data, index=[f'Size {size}' for size in sizes])
    #     csv_path = os.path.join(output_dir, 'time_comparison_data.csv')
    #     df.to_csv(csv_path)
    #     print(f"时间数据已保存到CSV文件: {csv_path}")

# 主函数
def main():
    # 读取多个数据集
    datasets = []
    
    # 假设有多个文件，分别替换为你的数据文件路径
    file_list = ["FISTA*", "PR_exp", "Elist++", "PR", "PR_lin", "Greedy++"]
    name = "_time.txt"
    filenames = []
    for each in file_list:
        filenames.append(each + name)
    
    for i, filename in enumerate(filenames):
        if os.path.exists(filename):
            sizes, averages = read_data(filename)
            color = file_colors.get(file_list[i], 'white')  # 获取颜色，如果没有则使用白色
            datasets.append((sizes, averages, file_list[i], color))
        else:
            print(f"警告: 文件 {filename} 不存在")
    
    # 将图表保存到 figures 目录下
    output_dir = 'figures'
    
    if datasets:
        plot_data(datasets, output_dir)
        
        # 打印数据摘要
        print("\n数据摘要:")
        for sizes, averages, name, color in datasets:
            print(f"{name}: 数据点数量={len(averages)}, 时间范围={min(averages):.4f}s - {max(averages):.4f}s")
    else:
        print("没有找到有效数据")

if __name__ == "__main__":
    main()
