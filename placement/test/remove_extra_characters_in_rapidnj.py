import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument('--dataset_size', type=int, help='dataset size')
args = parser.parse_args()

def remove_quotes_from_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()
    
    content_no_quotes = content.replace("'", "").replace('"', "")
    
    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(content_no_quotes)

file_path = f'/data/zec022/rapidNJ/dataset_{args.dataset_size}/tree.nwk'
remove_quotes_from_file(file_path)
