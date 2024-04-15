import numpy as np

print(str(np.where(np.array([0,1,1,0,0,1]))[0]))

print(", ".join(map(str, np.where(np.array([0,1,1,0,0,1]))[0]))+", ")
