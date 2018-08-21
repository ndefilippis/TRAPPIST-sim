from plot_closest import generate_mass_map

path = "distance_ecc"

data = generate_mass_map(path, contains_ecc=True)
for key in data.keys():
    for i in range(len(data[key])):
        data[key][i] = int(data[key][i][0].number) 
for ecc in data.keys():
    count = 0
    for mass in data[ecc]:
        for other_ecc in data.keys():
	    if other_ecc == ecc:
	        continue
            if mass not in data[other_ecc]:
	        count += 1
    print "Eccentricity", ecc, "has", count, "extra masses"
