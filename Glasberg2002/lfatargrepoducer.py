

class KV :

    def __init__(self):
        self.order = None
        self.fhigh = None
        self.flow = None
        order = None
        erbStep = None
        erbFcMin = None
        erbFcMax = None

    def setStandard(self):
            
        self.fs = 32000
        self.flow = 20
        self.fhigh = 16000
        self.order = 4096
        self.erbStep = 0.2500
        self.erbFcMin = 50
        self.erbFcMax = 15000