import os

path = "/home/draco/ndefilippis/distance_ecc"
directories = os.listdir(path)
list = ["-1", "1.0", "0.1", "0.01"]

for dir in directories:
#for n in range(0, 2000):
    prefix = dir[:9]
    #prefix = "close{0:03d}".format(n)
    for suffix in list:
        if prefix+suffix not in directories:
	    print prefix+suffix, "does not exist, while", dir, "does."
	#test_path = os.path.join(path, prefix+suffix)
	#if not os.path.exists(test_path):
	#    print test_path, "does not exist."

