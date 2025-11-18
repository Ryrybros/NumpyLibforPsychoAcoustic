
from jsonHandler import jsonHandler

import numpy as np


def FIR():
    
    data = jsonHandler.readJson("data/Glasberg2002.json")    
    tfLinear = 10.^(data["OuterMiddleEar1997"]["tfOuterMiddle"]/10)
    print(tfLinear)


FIR()

