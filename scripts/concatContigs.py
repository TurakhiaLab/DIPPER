from sys import argv
import os

path = "/data/swalia/drosophila_20/"
files=[]
for x in os.listdir(path):
    if (x.endswith(".fa")):
        files.append(x)

for f in files:
    print (f)
    filePath=path+f
    file=open(filePath, "r")
    data=file.readlines()
    file.close()

    filename=filePath[:-3] + "-concat.fa"
    file=open(filename, "w")

    file.write(">")
    file.write(f)
    file.write("\n")

    newData=""
    for line in data:
        if (line[0]==">"):
            continue
        else:
            file.write(line)

    file.close()

# for x in os.listdir(path):
#     if (x.endswith("concat.fa")):
        
