#code to find the potentially isolated points from all the data

data_files = ['1to100000.out', '100001to200000.out', '200001to300000.out', '300001to400000.out', '400001to500000.out', '500001to534299.out']

for name in data_files:
    file = open(name, 'r')
    data = file.readlines()
    bad_curves = []
    for line in data:
        jinv, ainv, b = line.split(':')
        if b.strip() != 'true':
            bad_curves.append(line)
    file.close()

if bad_curves:
    fileout = open('potentially_isolated.txt', 'w')
    for curve in bad_curves:
        fileout.write(curve)
    fileout.close()
else:
    print("No potentially isolated curves were found.")
