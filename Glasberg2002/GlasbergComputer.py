import numpy as np
from filetools import jsonHandler

Json = jsonHandler.jsonHandler
class GlasbergComputer:

    def __init__(self):
        self.data = Json.readJson("data/Glasberg2002.json")

    def OuterMiddle(self,fVec : np.array, model : str, free : bool ):

        # transfer function of the outer ear
        dataModel = f"OuterMiddleEar{model}"
        fOuter = np.array(self.data[dataModel]["fOuter"])
        

    #     % values of ANSI S3.4-2007
    #     if kv.fieldType == 'free' % free field
    #         tfOuter = [0 0 0 0 0 0 0 0 0.1 0.3 0.5 0.9 1.4 1.6 1.7 2.5 2.7 2.6 2.6 3.2 5.2 ...
    #             6.6 12 16.8 15.3 15.2 14.2 10.7 7.1 6.4 1.8 -0.9 -1.6 1.9 4.9 2 -2 2.5 2.5];

    #     elseif kv.fieldType == 'diffuse' % diffuse field
    #         tfOuter = [0 0 0 0 0 0 0 0 0.1 0.3 0.4 0.5 1 1.6 1.7 2.2 2.7 2.9 3.8 5.3 6.8 7.2 ...
    #             10.2 14.9 14.5 14.4 12.7 10.8 8.9 8.7 8.5 6.2 5 4.5 4 3.3 2.6 2 2];
    #     else
    #         error('Wrong parameter for fieldType. Please use "free" or "diffuse".')
    #     end

    #     tfOuterInterp = interp1(fOuter, tfOuter, kv.fVec, 'pchip');

    #     % transfer function of the middle ear
    #     fMiddle = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 750 ...
    #         800 1000 1250 1500 1600 2000 2500 3000 3150 4000 5000 6000 6300 8000 ...
    #         9000 10000 11200 12500 14000 15000 16000 18000 20000];

    #     % revised data 2006
    #     tfMiddle = -[39.6 32 25.85 21.4 18.5 15.9 14.1 12.4 11 9.6 8.3 7.4 6.2 4.8 ...
    #         3.8 3.3 2.9 2.6 2.6 3.2 4.5 5.5 8.5 10.4 7.3 7 6.6 7 9.2 10.2 12.2 ...
    #         10.8 10.1 12.7 15 18.2 23.8 32.3 45.5 50];

    #     tfMiddleInterp = interp1(fMiddle, tfMiddle, kv.fVec, 'pchip');

    #     data.tfOuterMiddle = tfOuterInterp + tfMiddleInterp;
    #     data.tfOuter = tfOuter;
    #     data.tfMiddle = tfMiddle;
    #     data.fOuter = fOuter;
    #     data.fMiddle = fMiddle;
    # end
    

    

if __name__ == '__main__':
    g = GlasbergComputer()
    g.OuterMiddle(",","2007",True) 