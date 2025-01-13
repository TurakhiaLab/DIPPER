import subprocess
import re
import random
import os
from ete3 import Tree


def replace_numbers(input_string):
    def replace_match(match):
        # 随机选择一个值
        return f":{random.choice([0.001, 0.0005, 0.0003, 0.0001, 0.0002, 0.0007, 0.0004, 0.0006, 0.0008, 0.0009])}"

    # 使用正则表达式匹配 ":" 后面的数字
    modified_string = re.sub(r':\d+(\.\d+)?', replace_match, input_string)
    return modified_string

def execute_command_realtime(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, universal_newlines=True)
    for line in iter(process.stdout.readline, ''):
        print(line, end='', flush=True)
    for line in iter(process.stderr.readline, ''):
        print(line, end='', flush=True)

def write_string_to_file(file_path, content):
    try:
        with open(file_path, 'w', encoding='utf-8') as file:
            file.write(content)
        print(f"Successfully wrote to the file: {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

def read_string_from_file(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()
        print(f"Successfully read from the file: {file_path}")
        return content
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def split_file(input_filename):
    with open(input_filename, 'r', encoding='utf-8') as infile:
        lines = infile.readlines()

    for i in range(0, len(lines), 2):
        if i + 1 < len(lines):
            title = lines[i].strip().replace(' ', '_').replace('/', '_')
            output_filename = f"./seqs/{title[1:]}.fa"
            with open(output_filename, 'w', encoding='utf-8') as outfile:
                outfile.write(lines[i])
                outfile.write(lines[i + 1])

def change_directory(new_directory):
    try:
        # 改变当前工作目录
        os.chdir(new_directory)
        print(f"Changed directory to: {new_directory}")
    except Exception as e:
        print(f"An error occurred: {e}")

def replace_first_line(filename, new_first_line):
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        if lines:
            lines[0] = new_first_line + '\n'
        else:
            lines.append(new_first_line + '\n')
        with open(filename, 'w', encoding='utf-8') as file:
            file.writelines(lines)
    
    except Exception as e:
        print(f"An error occurred: {e}")


def load_tree_from_file(filename):
    """
    Load a tree from a Newick format file.
    """
    try:
        with open(filename, 'r') as file:
            newick_str = file.read().strip()
            tree = Tree(newick_str)
            return tree
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


change_directory(os.path.expanduser('~/phylo-accel/test'))
t1 = load_tree_from_file("modified_tree.phy")
t2 = load_tree_from_file("phastsim_dataset/estimated_tree.nwk")
t3 = load_tree_from_file("phastsim_dataset/nj_tree.nwk")

# print(t1.robinson_foulds(t2))
result = t1.robinson_foulds(t2,unrooted_trees=True)
rf = result[0]
max_rf = result[1]
# print(t1, t2)
print("Our method:")
print("RF distance is %s over a total of %s" %(rf, max_rf))
print("Normalized RF is", rf/max_rf)
result = t1.robinson_foulds(t3,unrooted_trees=True)
rf = result[0]
max_rf = result[1]
# print(t1, t2)
print("Conventional Neighbor Joining method:")
print("RF distance is %s over a total of %s" %(rf, max_rf))
print("Normalized RF is", rf/max_rf)

# t4 = load_tree_from_file("phastsim_dataset/fastme_tree.nwk")
# result = t1.robinson_foulds(t4,unrooted_trees=True)
# rf = result[0]
# max_rf = result[1]
# # print(t1, t2)
# print("FastME method:")
# print("RF distance is %s over a total of %s" %(rf, max_rf))
# print("Normalized RF is", rf/max_rf)