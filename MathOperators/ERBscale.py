import numpy as np

class ERB :

    def f2erbrate(f : np.array, *erbModel : str):
        
        if(len(erbModel) < 1):
            model = 'moore1983'
        else:
            model = erbModel[0] 

        
        match model :
            case 'moore1983':
                f = f / 1000
                erbrate = 11.17 * np.log((f+0.312)/(f+14.675)) + 43

            case 'glasberg1990':
                f = f / 1000
                erbrate = 21.366 * np.log10(4.368*f + 1)
            case _:
                Exception('Unknown model for the conversion.')
                return

        
        return erbrate
        
    def erbrate2f(erbrate: np.array, *erbModel: str):
        #performs reverse transform as f2erbrate

        if(len(erbModel) < 1):
            model = 'moore1983'
        else:
            model = erbModel[0] 

        
        match model :
                case 'moore1983':
                    f = (0.312 - (np.exp((erbrate - 43)/11.17)) * 14.675) / (np.exp((erbrate - 43)/11.17) - 1)
                    f = f * 1000; 
                case 'glasberg1990':
                        f = (10**(erbrate/21.366)-1)/4.368
                        f = f * 1000
                case _:
                    Exception('Unknown model for the conversion.')
                    return
            
        return f
                
        


        
if __name__ == '__main__':
    
    max = 1000
    n = 100
    x = np.array([(max/n)*i for i in range(n)])
    try:
        scale = ERB.f2erbrate(x,'glasberg1990')
        f = ERB.erbrate2f(scale, 'glasberg1990')

    except :
        print("ERROR : The model is likely not correct")
        scale = x

    try:
        scale2 = ERB.f2erbrate(x,'moore1983')
    except :
        print("ERROR : The model is likely not correct")
        scale2 = x

    # print(scale)
    import matplotlib.pyplot as plt
    plt.plot(x,scale)
    plt.plot(x,f)
    # plt.plot(x,scale2)



    plt.show()