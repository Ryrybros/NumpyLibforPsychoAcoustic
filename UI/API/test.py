from Models.model import model
import numpy as np
with open("sine.txt", 'w') as f:
    y = [float(np.sin(5*i/32000) )for i in range(32000)]
    f.write(str(y))