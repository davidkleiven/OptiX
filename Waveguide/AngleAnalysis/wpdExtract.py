import numpy as np

class WpdDataset:
    def __init__( self ):
        self.x = None
        self.y = None
        self.name = ""

class WebPlotDigitizerDsetExtractor:
    def __init__( self ):
        self.dsets = []

    def parse( self, obj ):
        datasets = obj["wpd"]["dataSeries"]
        for dset in datasets:
            currentDset = WpdDataset()
            x = []
            y = []
            for val in dset["data"]:
                x.append( float(val["value"][0]) )
                y.append( float(val["value"][1]) )
            currentDset.x = np.array( x )
            currentDset.y = np.array( y )
            currentDset.name = dset["name"]
            self.dsets.append(currentDset)

    def get( self, name ):
        for dset in self.dsets:
            if ( dset.name == name ):
                return dset
        print ("Warning! Did not find any dataset with name: %s"%(name))
        return None
