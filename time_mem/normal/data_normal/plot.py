import numpy as np
import matplotlib.pyplot as plt
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

# 画图并保存为 PDF 和 PNG 文件
def plot_data(datasets, output_dir):
    # 创建 figures 目录，如果不存在的话
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 图像文件路径
    pdf_path = os.path.join(output_dir, 'output.pdf')
    png_path = os.path.join(output_dir, 'output.png')
    eps_path = os.path.join(output_dir, 'output.eps')
    plt.figure(figsize=(10, 6))  # 设置图形大小
    # 设置 x 轴和 y 轴为对数刻度
    plt.xscale('log')
    plt.yscale('log')
    #colors = ['blue', 'green', 'orange', 'purple', 'cyan', 'magenta']  # 可用的颜色列表
    # 画不同数据集的折线图
    for i, (sizes, averages, label, color) in enumerate(datasets):
        #color = colors[i % len(colors)]  # 为每个数据集分配颜色
        plt.plot(sizes, averages, linestyle='-', marker='o', color=color, label=label)

    # 设置图表标题和轴标签
    #plt.title('Comparison of Average Time vs Size for Different Data Sets')
    plt.xlabel(r'Size $|\mathcal{F}|$')  # 使用 LaTeX 格式
    plt.ylabel('Average Time per iteration (second)')

    # 显示图例
    plt.legend()

    # 保存图表为 PDF 文件
    plt.savefig(pdf_path)
    print(f'图表已保存为 {pdf_path}')

    # 保存图表为 PNG 文件
    plt.savefig(png_path)
    print(f'图表已保存为 {png_path}')

    # 保存图表为 eps 文件
    plt.savefig(eps_path)
    print(f'图表已保存为 {eps_path}')
    # 关闭图表显示，防止内存泄漏
    plt.close()


# 定义文件的颜色映射
file_colors = {
    'PR_lin': 'magenta',#'pink'
    'PR': 'orange',
    'PR_exp': 'red',
    # 你可以继续为其他文件手动设置颜色
    #'FISTA': 'blue',
    'FW_Elist': 'lightblue',
    'Elist++': 'green',
    'Greedy++': 'darkblue',
    #'FISTA': 'cyan',
    'FISTA*': 'cyan'
    # 其他文件可以继续添加...
}

# 主函数
def main():
    # 读取多个数据集
    datasets = []
    
    # 假设有多个文件，分别替换为你的数据文件路径
    file_list = ["FISTA*", "PR_exp","Elist++", "PR", "PR_lin", "Greedy++"] # 
    name = "_time.txt"
    filenames = []
    for each in file_list:
        filenames.append(each+name)
    
    for i, filename in enumerate(filenames):
        sizes, averages = read_data(filename)
        datasets.append((sizes, averages, file_list[i], file_colors[file_list[i]]))

    # 将图表保存到 figures 目录下
    output_dir = 'figures'
    plot_data(datasets, output_dir)

main()
