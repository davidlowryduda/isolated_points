for x in ["1to100000","100001to200000", "200001to300000", "300001to400000", "500001to534299", "genusgt0"]:  
	with open(x) as f:
		lines = f.readlines()
		for l in lines:
			j, pts = l.split(":")
			if len(pts) > 6:
				print(j)
				print(pts)
				