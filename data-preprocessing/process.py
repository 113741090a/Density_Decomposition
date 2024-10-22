# -*- coding: utf-8 -*-


def print_lines(file_path, x, y):
    with open(file_path, 'r') as file:
        # ��ӡ�ļ���ǰ x ��
        for _ in range(x):
            print(file.readline().rstrip())

        # ��ȡ�ļ�������������λ��ĩβ y ��֮ǰ��λ��
        file.seek(0, 2)
        end_pos = file.tell()
        file.seek(0, 0)
        total_lines = sum(1 for _ in file)

        start_pos = 0
        if total_lines > y:
            for _ in range(total_lines - y):
                start_pos = file.readline()

        # ��ӡ�ļ��ĺ� y ��
        file.seek(start_pos)
        for line in file:
            print(line.rstrip())


def delete_lines(file_path, x):
    # ɾ���ļ���ǰ x ��
    with open(file_path, 'r+') as file:
        lines = file.readlines()
        del lines[:x]

    with open(file_path, 'w') as file:
        file.writelines(lines)

from collections import deque

def print_first_and_last_lines(filename, x, y):
    with open(filename, 'r') as file:
        # ��ȡǰ x ��
        print(f"First {x} lines:")
        for i in range(x):
            line = file.readline()
            if not line:
                break  # �ļ���������x��
            print(line.strip())

        print(f'\nLast {y} lines:')
        # ��ȡ�ļ������ y ��
        # deque ֻ������� y ������
        last_lines = deque(file, maxlen=y)
        for line in last_lines:
            print(line.strip())

def delete_first_x_lines(filename, x):
    with open(filename, 'r') as file:
        # ����ǰ x ��
        for _ in range(x):
            next(file, None)
        
        # ��ʣ����б���
        remaining_lines = file.readlines()
    
    # д��ʣ����У�����ԭ�ļ�
    with open(filename, 'w') as file:
        file.writelines(remaining_lines)

# ���ú�����ָ���ļ�·����Ҫ��ӡ��ǰ x �кͺ� y ��
file_path = 'graph.txt'
print_first_and_last_lines(file_path, 5, 3)
#delete_first_x_lines(file_path, 4)
#print_first_and_last_lines(file_path, 5, 3)
