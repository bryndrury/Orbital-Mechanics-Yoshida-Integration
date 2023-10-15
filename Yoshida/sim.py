import os
import time
from animate import animate


os.system("clang++ -std=c++17 -Ofast -march=native vector3d.cpp body.cpp compute_orbits.cpp -o orbit_sim")
#os.system("clang++ -std=c++17 -g vector3d.cpp body.cpp compute_orbits.cpp -o orbit_sim")

start = time.time()
# 40750 * 2
os.system("./orbit_sim solar_system.csv solar_results.csv 1e-6 815000 750")
print(f"Time taken: {time.time()-start}s")

# os.system("python3 postprocess.py")
# print("Post Processing Completed.")

# animate("solar_results", 30)
# animate("solar_results", 5)
