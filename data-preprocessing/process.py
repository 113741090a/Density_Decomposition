# -*- coding: utf-8 -*-
import sys




def print_lines(file_path, x, y):
    with open(file_path, 'r') as file:
        for _ in range(x):
            print(file.readline().rstrip())

        file.seek(0, 2)
        end_pos = file.tell()
        file.seek(0, 0)
        total_lines = sum(1 for _ in file)

        start_pos = 0
        if total_lines > y:
            for _ in range(total_lines - y):
                start_pos = file.readline()


        file.seek(start_pos)
        for line in file:
            print(line.rstrip())


def delete_lines(file_path, x):

    with open(file_path, 'r+') as file:
        lines = file.readlines()
        del lines[:x]

    with open(file_path, 'w') as file:
        file.writelines(lines)

from collections import deque

def print_first_and_last_lines(filename, x, y):
    with open(filename, 'r') as file:

        print(f"First {x} lines:")
        for i in range(x):
            line = file.readline()
            if not line:
                break  
            print(line.strip())

        print(f'\nLast {y} lines:')

        last_lines = deque(file, maxlen=y)
        for line in last_lines:
            print(line.strip())

def delete_first_x_lines(filename, x):
    with open(filename, 'r') as file:

        for _ in range(x):
            next(file, None)
        

        remaining_lines = file.readlines()
    

    with open(filename, 'w') as file:
        file.writelines(remaining_lines)

file_path = 'graph.txt'

# sys.argv[0] 是脚本名称
# sys.argv[1] 是第一个参数，以此类推
if len(sys.argv) > 1:
    print("接收到的参数：", sys.argv[1:])
    num_lines = int(sys.argv[1])
    delete_first_x_lines(file_path, 4)
else:
    print("没有接收到参数")
    #delete_first_x_lines(file_path, 4)

#print_first_and_last_lines(file_path, 5, 3)
#delete_first_x_lines(file_path, 4)
#print_first_and_last_lines(file_path, 5, 3)
