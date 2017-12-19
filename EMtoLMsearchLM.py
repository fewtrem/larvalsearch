'''
Created on 11 Nov 2017

@author: fewtrem
'''
# let's calculate the mins first!
from constants import BRAIN
import numpy as np,os

import ctypes
curDir = os.path.abspath(os.path.dirname(__file__))
lib = ctypes.cdll.LoadLibrary(os.path.join(curDir,"compareSkeletonPointsProd2.so"))
getNBlastScore = lib.getNBlastScoreProd
# get the minimum distance for all points in C1 to those C2 lists:
def cGetNBlastScore(testC1,testC2,compDataA,compDataB):
    #Order C2 to search for points
    cGetOrderIndecies = lib.getOrderIndecies
    testB = np.ascontiguousarray(testC2.astype(np.intc))
    testBShape = np.array(testB.shape).astype(np.intc)
    outputB = np.zeros((1000+1))
    outputB = np.ascontiguousarray(outputB.astype(np.intc))
    outputBShape = np.array(outputB.shape).astype(np.intc)
    cGetOrderIndecies(ctypes.c_void_p(testB.ctypes.data),
             ctypes.c_void_p(testBShape.ctypes.data),
             ctypes.c_void_p(outputB.ctypes.data),
             ctypes.c_void_p(outputBShape.ctypes.data))
    #Next
    testCO = outputB
    testCO = np.ascontiguousarray(testCO.astype(np.intc))
    testC1 = np.ascontiguousarray(testC1.astype(np.intc))
    testC1Shape = np.array(testC1.shape).astype(np.intc)
    testC2 = np.ascontiguousarray(testC2.astype(np.intc))
    testC2Shape = np.array(testC2.shape).astype(np.intc)
    compDataA = np.ascontiguousarray(compDataA.astype(np.float32))
    compDataB = np.ascontiguousarray(compDataB.astype(np.float32))
    outputC = np.zeros((testC1.shape[0],4),dtype=np.float32)
    outputC = np.ascontiguousarray(outputC.astype(np.float32))
    getNBlastScore(ctypes.c_void_p(testCO.ctypes.data),
             ctypes.c_void_p(testC1.ctypes.data),
             ctypes.c_void_p(testC1Shape.ctypes.data),
             ctypes.c_void_p(compDataA.ctypes.data),
             ctypes.c_void_p(testC2.ctypes.data),
             ctypes.c_void_p(testC2Shape.ctypes.data),
             ctypes.c_void_p(compDataB.ctypes.data),
             ctypes.c_void_p(outputC.ctypes.data))
    # Return is DistanceSquared,Dot Product,PositionSquared BLAST Result, NBlast Result.
    return  outputC


def getAllForOne(mU150,testC1,compDataA):
    thisType = BRAIN
    # reset this value's 1000s to something else so they don't equate!
    testC1[testC1==10000]=30000
    selectorAB = testC1[:,1]<10000
    thisResList = []
    for jKey in range(mU150.pointsC[thisType].shape[0]):
        if jKey%10000==0:
            print jKey
        #if goodOrNot[jKey]==True:
        testC2 = mU150.pointsC[thisType][jKey]
        selectorBA = testC2[:,1]<10000
        if np.sum(selectorBA)>0:
            compDataB = mU150.nBlastC[thisType][jKey]
            thisRes = cGetNBlastScore(testC1,testC2,compDataA,compDataB)
            compResAB = thisRes[selectorAB]
            thisResInv = cGetNBlastScore(testC2,testC1,compDataB,compDataA)
            compResBA = thisResInv[selectorBA]
            thisPredABNB = np.sum(compResAB[:,3])*1.0/np.sum(selectorAB)
            thisPredBANB = np.sum(compResBA[:,3])*1.0/np.sum(selectorBA)
            thisPredABPB = np.sum(compResAB[:,2])*1.0/np.sum(selectorAB)
            thisPredBAPB = np.sum(compResBA[:,2])*1.0/np.sum(selectorBA)
            thisResList.append([jKey,thisPredABNB,thisPredBANB,thisPredABPB,thisPredBAPB])
    return thisResList