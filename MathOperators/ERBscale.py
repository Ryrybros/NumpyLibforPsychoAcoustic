import numpy as np

class ERB :

    # This class has all static methods that are used in the model that are realted to ERB scales, the names match the amt-toolbox names

    def f2erbrate(f : np.array, erbModel = "glasberg1990"):
        
        model = erbModel
        
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
        
    def erbrate2f(erbrate: np.array, erbModel = 'glasberg1990' ):
        #performs reverse transform as f2erbrate
        model = erbModel
        
        match model :
                case 'moore1983':
                    f = (0.312 - (np.exp((np.array(erbrate) - 43)/11.17)) * 14.675) / (np.exp((np.array(erbrate) - 43)/11.17) - 1)
                    # f = (0.312 - (np.exp((erbrate - 43*np.ones(len(erbrate)) )/11.17)) * 14.675) / (np.exp((erbrate - 43*np.ones(len(erbrate)) )/11.17) - 1)
                    f = f * 1000; 
                case 'glasberg1990':
                        f = (10**(erbrate/21.366)-1)/4.368
                        f = f * 1000
                case _:
                    Exception('Unknown model for the conversion.')
                    return
        # print("f of erbrate2f is : ", f)
        return f
    

        

    
    def f2erb( f : float,
               ERB_different_values  = False
                   ):

        if( ERB_different_values== False  ):
            return 24.673*(0.004368*f + 1) 
        else:
            ERB_Q = 1000/(24.7*4.37) # 9.2645
            
            ERB_f = 1000/4.37 # 228.833

        return (ERB_f + f) / ERB_Q
        
        
        
        


        
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

    
    import matplotlib.pyplot as plt
    plt.plot(x,scale)
    
    # plt.plot(x,f)

    plt.plot(x,scale2)



    plt.show()

    print(ERB.f2erb(0))