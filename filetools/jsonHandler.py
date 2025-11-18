
import json

class jsonHandler:
            
    def readJson(path : str):
        with open(path, 'r') as file:
            python_obj = json.load(file)
        return python_obj
    
    def writeJson(filePath : str, data : dict):
        
        try:
            with open(filePath, 'w') as file:
                json.dump(data, file)
            print("Succesfully saved json file")
        except:
            print("Failed to save json file")

#test
if __name__ == '__main__':
    data = jsonHandler.readJson("filetools/data.json")
    assert(data["OuterMiddleEar1997"]["fOuter"][0] == 20)
    jsonHandler.writeJson("filetools/testjsonHandler.json", data)
            
        