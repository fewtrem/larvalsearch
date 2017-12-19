'''
Created on 8 Dec 2017

SNAPSHOT VERSION WITH DIRECTORIES CHANGED!

@author: fewtrem

# run using:
java -jar jython.jar "/afs/inf.ed.ac.uk/user/s11/s1144899/PhD/Python Projects/TrackEMToSVG/Eclipse/runLandmarkReg.py"


'''
import sys,os
import cPickle as pickle
curDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append("/plugin/mpicbg-master/mpicbg/target/mpicbg-1.1.2-SNAPSHOT.jar")
from mpicbg.models import Point, PointMatch, AffineModel3D, MovingLeastSquaresTransform

with open(os.path.join(curDir,"MLS","pointsToMatchDec2017EntryPointsB.pkl")) as fI:
    landmarkList = pickle.load(fI) 
matches = []
for thisLand in landmarkList:
    for thisSide in ["right","left","centre"]:
        if "FROM_"+thisSide in thisLand and "TO_"+thisSide in thisLand:
            p1 = Point(map(float,thisLand["FROM_"+thisSide]))
            p2 = Point(map(float,thisLand["TO_"+thisSide]))
            matches.append(PointMatch(p1, p2, 1.0))

print len(matches),"landmark points used"

mls = MovingLeastSquaresTransform()
model = AffineModel3D()
mls.setModel(model)
mls.setMatches(matches)

with open(os.path.join(curDir,"MLS","Points_input.pkl")) as fI:
    pointList = pickle.load(fI)
output = []
for i in range(len(pointList)):
        output.append(map(float, mls.apply(pointList[i])))

with open(os.path.join(curDir,"MLS","ConvertedPoints_output.pkl"),'w') as fO:
    pickle.dump(output,fO)