import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument('--inp', type=str, help='input file path')
parser.add_argument('--out', type=str, help='out file path')
parser.add_argument('--type', type=str, help='quicktree/rapidnj')
args = parser.parse_args()

def remove_newline_from_file(file_path, output_file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()
    
    content_no_newline_dot = content.replace('\n', '')
    
    with open(output_file_path, 'w', encoding='utf-8') as file:
        file.write(content_no_newline_dot)

def remove_quotes_from_file(file_path, output_file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()
    
    content_no_quotes = content.replace("'", "").replace('"', "")
    
    with open(output_file_path, 'w', encoding='utf-8') as file:
        file.write(content_no_quotes)

if (args.type == "quicktree"):
    remove_newline_from_file(args.inp, args.out)
elif (args.type == "rapidnj"):
    remove_quotes_from_file(args.inp, args.out)