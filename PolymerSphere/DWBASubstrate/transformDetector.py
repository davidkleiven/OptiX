import numpy as np

class DetectorTransformer:
    def __init__(self):
        self.name = "identity"
    
    def scatteringAngle(self, incidentGrazingAngle, grazingAngleFromSubstrate):
        return grazingAngleFromSubstrate

    def axisLabel(self):
        return "$\\alpha_f/\\alpha_c$"

class DetectorCenterBeam(DetectorTransformer):
    def __init__(self):
        self.name = "detectorCenterBeam"

    def scatteringAngle(self, incidentGrazingAngle, grazingAngleFromSubstrate):
        return incidentGrazingAngle + grazingAngleFromSubstrate
    
    def axisLabel(self):
        return "$(\\alpha+\\alpha_f)/\\alpha_c$"
         
