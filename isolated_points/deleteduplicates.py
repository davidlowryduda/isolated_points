#code to get isolated_refined.txt from isolated.txt

file = open('isolated.txt', 'r')
data = file.readlines()
from collections import defaultdict
dict = defaultdict(list)
for line in data:
    jinv, ainv = line.split('","')
    a,b =jinv[1:].split(',')
    j = (eval(a[1:]),eval(b[:-1]))
    ainv = ainv[:-2].replace('{','[')
    ainv = ainv.replace('}', ']')
    ainvs = eval(ainv)
    dict[j].append(ainv)
file.close()

fileout = open('isolated_refined.txt', 'w')
for j in dict.keys():
    fileout.write(str(j))
    fileout.write(',')
    fileout.write(str(dict[j][0]))
    fileout.write('\n')
fileout.close()