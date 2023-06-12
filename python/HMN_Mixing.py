from HiggsAnalysis.CombinedLimit.PhysicsModel import *

class HMN_Mixing_Steering(PhysicsModel):
    def __init__(self):
        self.pois = {}
    
    def setModelBuilder(self,modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False
    
    def getYieldScale(self,bin,process):
        if 'hmn_m' in process:
            return "r * r"
        elif:
            'HN' in process:
            return "r"
        else:
            return 1
    def doParametersOfInterest(self):
        self.modelBuilder.doVar("r[0.,-100.,100.]"):
        self.modelBuilder.doSet("POI","r")

hmn_mixing_steering = HMN_Mixing_Steering
