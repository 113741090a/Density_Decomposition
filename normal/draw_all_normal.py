import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import FuncFormatter, FixedLocator
#algo_list = [ "FISTA", "FW_Elist", "PR", "PR_lin", "Greedypp", "Elistpp","PR_exp_pro"]
algo_names = ["FISTA*", "Greedy++", "Elist++", "PR", "PR_lin", "PR_exp"]

# 自定义y轴刻度的函数
def custom_scale(y, pos):
    if y <= 10:
        return f'{y:.0f}'  # 0-10 线性刻度
    else:
        # 10以上采用对数刻度标签
        return f'$10^{{{int(np.log10(y))}}}$'

def scientific_notation(num):
    if num == 0:
        return (0, 0)
    
    exponent = 0
    while abs(num) >= 10:
        num /= 10
        exponent += 1
    
    while abs(num) < 1:
        num *= 10
        exponent -= 1

    return (num, exponent)

def calculate_l2_norm(file_path):
    res = 0.0
    coordinates = []
    with open(file_path, 'r') as file:
        for line in file:
            numbers = line.strip().split()
            if len(numbers) >= 2:
                res = res + float(numbers[1])*float(numbers[1])

    
    res = res**(1/2)

    return res

file_path = 'Exact_normal.txt'  # 替换为你的文件路径

real_obj = calculate_l2_norm(file_path)
res = scientific_notation(real_obj)
num = res[0]
exponent = res[1]

# 自定义刻度
# 自定义的函数，用于将原本的10^y替换为10^(y-x)
def custom_format(y, pos, x=3):
    """将10^y转换成10^(y-x)的形式"""
    exponent = int(np.log10(y))  # 获取log scale的指数部分
    return r'$10^{{{}}}$'.format(exponent - x)

def divide_elements(arr, constant):
    new_arr = [elem / constant for elem in arr]
    return new_arr

def read_data(file, swap_columns=False, skip=0):
    """
    读取文件中的数据，当 `swap_columns` 为 True 时，交换列顺序。
    """
    iterations = []
    errors = []
    with open(file, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count <= skip:
                continue
            data = line.split()
            if swap_columns:
                error, iteration = data
            else:
                iteration, error = data
            iterations.append(float(iteration))
            errors.append(float(error))
    return iterations, errors


def plot_data(files, title, xlabel, ylabel, output_file, swap_columns=False, skip=0, log_scale=True, markersize=4, colors=None, normalized = False):
    """
    遍历文件，读取数据并绘制图形。
    
    - `markersize` 控制点的大小。
    - `colors` 是一个字典，key 是文件名，value 是对应的颜色。
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    count = 0
    for file in files:
        iterations, errors = read_data(file, swap_columns=swap_columns, skip=skip)
        #color = colors.get(file, None) if colors else None  # 获取对应的颜色，若未指定则使用默认颜色
        color = None
        length = 0
        for key, value in colors.items():
            if key in file:
                if (len(key) >= length):
                    color = value
                    length = len(key)
        if normalized:
            errors = divide_elements(errors, num)
                
        ax.plot(iterations, errors, marker='o', linestyle='-', label=algo_names[count], markersize=markersize, color=color)
        count = count + 1
    if log_scale:
        ax.set_yscale('log')
        if normalized:
            # 动态设置x的值
            x_dynamic = exponent  # 你可以根据需要动态调整这个值

            # 使用lambda函数传递动态的x参数
            ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: custom_format(y, pos, x_dynamic)))
            #ax.yaxis.set_major_formatter(FuncFormatter(custom_format))
            
        #plt.yscale('log')
    else:
        ax.set_yscale('symlog', linthresh=10)  # linthresh=10表示在10以下的部分为线性刻度
        ax.yaxis.set_major_formatter(FuncFormatter(custom_scale))

    # 添加标题和标签
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_title(title)
    # ax.title('Inversion Number vs Iterations')
    # ax.xlabel('Iterations')
    # ax.ylabel('Inversion Number')

    # 显示网格
    ax.grid(True, which="both", ls="--")

    # 显示图例
    ax.legend()        
    # plt.title(title)
    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    # plt.grid(True, which="both", ls="--")
    # plt.legend()

    # 保存为 PNG 文件
    plt.savefig(output_file + '.png', dpi=300)
    # 保存为 PDF 文件
    pdf_path = "pdf/"
    pdf_path = ""
    plt.savefig(pdf_path+output_file + '.pdf')  # PDF 通常不需要 dpi 设置，因为是矢量图
    plt.savefig(pdf_path+output_file + '.eps', format = "eps")
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
    'FISTA': 'cyan'
    # 其他文件可以继续添加...
}

# 定义绘图任务
folder_name = "data_normal/"
algo_list = ["FISTA", "Greedy++", "Elist++", "PR", "PR_lin", "PR_exp"]
#algo_list = ["FISTA", "Greedy++", "Elist++", "PR", "PR_lin", "PR_exp"]
#algo_list = [ "FISTA", "FW_Elist", "PR", "PR_lin", "Greedypp", "Elistpp","PR_exp_pro"] "FISTA*",
#"FW_syn","FW_WF","PR_exp",
mul_name = "_mul.txt"
mul_time_name = "_mul_time.txt"
inv_name = "_inv.txt"
inv_time_name = "_inv_time.txt"
abs_name = "_abs.txt"
abs_time_name = "_abs_time.txt"
abs_normalized_name = "_abs_normalized.txt"
abs_normalized_time_name = "_abs_normalized_time.txt"
obj_name = "_obj.txt"
obj_time_name = "_obj_time.txt"
obj_norm_name = "_obj_normalized.txt"
obj_norm_time_name = "_obj_normalized_time.txt"
#mul_time_name = "_mul_time.txt"

tasks = [
    {
        'files': [folder_name+algo + abs_name for algo in algo_list],
        'title': 'Global Error vs Iterations',
        'xlabel': 'Iterations',
        'ylabel': 'Global Error',
        'output_file': 'figures/Absolute_Error_vs_T',
        'swap_columns': False,
        'colors': file_colors,
        'normalized': True
    },
    {
        'files': [folder_name+algo + abs_time_name for algo in algo_list],
        'title': 'Global Error vs Time',
        'xlabel': 'Time (seconds)',
        'ylabel': 'Global Error',
        'output_file': 'figures/Absolute_Error_vs_Time',
        'swap_columns': True,
        'colors': file_colors,
        'normalized': True
    },
    {
        'files': [folder_name+algo + abs_normalized_name for algo in algo_list],
        'title': 'Absolute Error vs Iterations',
        'xlabel': 'Iterations',
        'ylabel': 'Error',
        'output_file': 'figures/Normalized_Absolute_Error_vs_T',
        'swap_columns': False,
        'colors': file_colors,
        'normalized': False
    },
    {
        'files': [folder_name+algo + inv_name for algo in algo_list],
        'title': 'Number of Inversions vs Iterations',
        'xlabel': 'Iterations',
        'ylabel': 'Number of Inversions',
        'output_file': 'figures/inv_vs_T',
        'swap_columns': False,
        'log_scale': False,
        'colors': file_colors,
        'normalized': False
    },
    {
        'files': [folder_name+algo + inv_time_name for algo in algo_list],
        'title': 'Number of Inversions vs Time',
        'xlabel': 'Time (seconds)',
        'ylabel': 'Number of Inversions',
        'output_file': 'figures/inv_vs_Time',
        'swap_columns': True,
        #'swap_columns': False,
        'log_scale': False,
        'colors': file_colors,
        'normalized': False
    },
    {
        'files': [folder_name+algo + mul_name for algo in algo_list],
        'title': 'Local Error vs Number of Iterations',
        'xlabel': 'Iterations',
        'ylabel': 'Local Error',
        'output_file': 'figures/Multiplicative_Error_vs_T',
        'swap_columns': False,
        'colors': file_colors,
        'normalized': False
    },
    {
        'files': [folder_name+algo + mul_time_name for algo in algo_list],
        'title': 'Local Error vs Time',
        'xlabel': 'Time (seconds)',
        'ylabel': 'Local Error',
        'output_file': 'figures/Multiplicative_Error_vs_Time',
        'swap_columns': True,
        'colors': file_colors,
        'normalized': False
    }
]

# 执行绘图
for task in tasks:
    plot_data(
        files=task['files'],
        title=task['title'],
        xlabel=task['xlabel'],
        ylabel=task['ylabel'],
        output_file=task['output_file'],
        swap_columns=task.get('swap_columns', False),
        skip=task.get('skip', 0),
        log_scale=task.get('log_scale', True),
        markersize=1,  # 设置较小的点大小
        colors=task.get('colors', None),
        normalized=task.get('normalized', True),
    )

print("Finished")


'''
    {
        'files': [folder_name+algo + abs_normalized_time_name for algo in algo_list],
        'title': 'Absolute Error vs Time',
        'xlabel': 'Time',
        'ylabel': 'Error',
        'output_file': 'figures/Normalized_Absolute_Error_vs_Time',
        'swap_columns': True,
        'colors': file_colors,
        'normalized': False
    },
    {
        'files': [folder_name+algo + inv_name for algo in algo_list],
        'title': 'Number of Inversions vs Iterations',
        'xlabel': 'Iterations',
        'ylabel': 'Number of Inversions',
        'output_file': 'figures/inv_vs_T',
        'swap_columns': False,
        'log_scale': False,
        'colors': file_colors,
        'normalized': False
    },
    {
        'files': [folder_name+algo + abs_normalized_name for algo in algo_list],
        'title': 'Absolute Error vs Iterations',
        'xlabel': 'Iterations',
        'ylabel': 'Error',
        'output_file': 'figures/Absolute_Error_vs_T',
        'swap_columns': False,
        'colors': file_colors
    },
    {
        'files': [folder_name+algo + abs_normalized_time_name for algo in algo_list],
        'title': 'Absolute Error vs Time',
        'xlabel': 'Time',
        'ylabel': 'Error',
        'output_file': 'figures/Absolute_Error_vs_Time',
        'swap_columns': True,
        'colors': file_colors
    },
def plot_data(files, title, xlabel, ylabel, output_file, swap_columns=False, skip=0, log_scale=True, markersize=4, colors=None, normalized = False):
    """
    遍历文件，读取数据并绘制图形。
    
    - `markersize` 控制点的大小。
    - `colors` 是一个字典，key 是文件名，value 是对应的颜色。
    """
    if (log_scale):
        fig, ax = plt.subplots()
        count = 0
        for file in files:
            iterations, errors = read_data(file, swap_columns=swap_columns, skip=skip)
            #color = colors.get(file, None) if colors else None  # 获取对应的颜色，若未指定则使用默认颜色
            color = None
            length = 0
            for key, value in colors.items():
                if key in file:
                    if (len(key) >= length):
                        color = value
                        length = len(key)
                    

            ax.plot(iterations, errors, marker='o', linestyle='-', label=algo_list[count], markersize=markersize, color=color)
            count = count + 1
        if log_scale:
            ax.set_yscale('log')
            #plt.yscale('log')
        # 添加标题和标签
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        # ax.title('Inversion Number vs Iterations')
        # ax.xlabel('Iterations')
        # ax.ylabel('Inversion Number')

        # 显示网格
        ax.grid(True, which="both", ls="--")

        # 显示图例
        ax.legend()        
        # plt.title(title)
        # plt.xlabel(xlabel)
        # plt.ylabel(ylabel)
        # plt.grid(True, which="both", ls="--")
        # plt.legend()
    else:
        fig, ax = plt.subplots()
        count = 0
        # 遍历每个文件，读取数据并绘制图形
        for file in files:
            iterations, errors = read_data(file, swap_columns=swap_columns, skip=skip)
            color = None
            length = 0
            for key, value in colors.items():
                if key in file:
                    if (len(key) >= length):
                        color = value
                        length = len(key)

            # 绘制每个文件的折线图
            ax.plot(iterations, errors, marker='o', linestyle='-',label=algo_list[count], markersize=markersize, color=color)
            count = count + 1
            #plt.plot(errors, iterations, marker='o', linestyle='-', label=file)

        # 设置对数刻度
        #plt.xscale('log')
        # 设置自定义的y轴刻度
        ax.set_yscale('symlog', linthresh=10)  # linthresh=10表示在10以下的部分为线性刻度
        ax.yaxis.set_major_formatter(FuncFormatter(custom_scale))


        # 添加标题和标签
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        # ax.title('Inversion Number vs Iterations')
        # ax.xlabel('Iterations')
        # ax.ylabel('Inversion Number')

        # 显示网格
        ax.grid(True, which="both", ls="--")

        # 显示图例
        ax.legend()

    # 保存为 PNG 文件
    plt.savefig(output_file + '.png', dpi=300)
    # 保存为 PDF 文件
    pdf_path = "pdf/"
    pdf_path = ""
    plt.savefig(pdf_path+output_file + '.pdf')  # PDF 通常不需要 dpi 设置，因为是矢量图

    plt.close()
    '''
