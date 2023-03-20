# code to get isolated_refined.txt from isolated1.txt and isolated2
from collections import defaultdict


with open('isolated1.txt', 'r', encoding="utf8") as file1:
    data1 = file1.readlines()
with open('isolated2.txt', 'r', encoding="utf8") as file2:
    data2 = file2.readlines()

data = data1 + data2

# remove comments at start of file
data = [datum for datum in data if datum[0] != "#"]
ndata = len(data)


print("Parsing isolated1 and isolated2 files...")
dict = defaultdict(list)
for idx, line in enumerate(data):
    jinv, ainv = line.split('","')
    a, b = jinv[1:].split(',')
    j = (eval(a[1:]), eval(b[:-1]))
    ainv = ainv[:-2].replace('{', '[')
    ainv = ainv.replace('}', ']')
    ainvs = eval(ainv)
    dict[j].append(ainv)
    print(f"\r{idx:10}/{ndata}", end="")
print("Done.")


with open("isolated_refined.txt", 'w', encoding="utf8") as fileout:
    for j in dict.keys():
        fileout.write(str(j))
        fileout.write(',')
        fileout.write(str(dict[j][0]))
        fileout.write('\n')
print("Refined output written to isolated_refined.txt")
