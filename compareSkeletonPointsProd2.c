//gcc -shared -fPIC -I/usr/include/python2.7/ -lpython2.7 -o compareSkeletonPointsProd2.so compareSkeletonPointsProd2.c
#include <math.h>


int doSum(const void * condArrP, const void * outDDimsP) {
    // Process the inputs:
    int * condArr = (int *) condArrP;
    int * outDDims = (int *)outDDimsP;
    // set up the variables we will use:
	int x,y,z,posX,posY,pos;
	int t;
	t=0;
	// go through the values in the outData:
	for(x=0;x<outDDims[0];x++){
		posX = x*outDDims[2]*outDDims[1];
		for(y=0;y<outDDims[1];y++){
			posY = y*outDDims[2];
			for (z=0;z<outDDims[2];z++){
				pos = posX+posY+z;// Position we will save in the outD.
				t+= condArr[pos];
			}
		}
	}
	return t;


}

int getOrderIndecies(const void * inputArrP, const void * inputArrDimsP,const void * outputBP,const void * outputBDimsP) {
    // Process the inputs:
    int * inputArr = (int *) inputArrP;
    int * inputArrDims = (int *) inputArrDimsP;
    int * outputB = (int *) outputBP;
    int * outputBDims = (int *) outputBDimsP;
    // set up the variables we will use:
	int x,posX,j;
	int lastOne,thisOne;
	lastOne=-1;
	// go through the values in the outData:
	for(x=0;x<inputArrDims[0];x++){
		posX = x*inputArrDims[1];
		thisOne = inputArr[posX];
		if (lastOne!=thisOne){
			for (j=lastOne+1;j<thisOne+1;j++){
				outputB[j]=x;
			}
		}
		
		lastOne = thisOne;
	}
	
	for (j = lastOne+1;j<outputBDims[0];j++){
		outputB[j]=inputArrDims[0];
	}
	return 0;
}



int getNBlastScoreProd(const void * inputArrCOP,
const void * inputArrAP,const void * inputArrADimsP,const void * condArrAP,
const void * inputArrBP, const void * inputArrBDimsP, const void * condArrBP,
const void * outputCP) {
	// constants:
	const float posThreshold = 25;
	const float posScoreIfOverThreshold = -3.25;
	const float nBlastThreshold= posThreshold;
	const float nBlastScoreIfOverThreshold = -3.25;
	const int threshhold = 87;
	const int sizeOfComp = 3;
	const int xAxis = 1000;
	float mult[3] = {0.293,0.293,0.5};
	// c+x+y+x**2+x*y+y**2+x**3+x**2y+xy**2+y**3, x dist squared, y dot product
	float nBCo[10] = {1.19641997e+00,-3.40663721e-01,2.60766902e+00,1.10076885e-02,-3.64288080e-01,8.19858762e-01,-1.73581331e-04,4.08108191e-03,2.27248084e-01,-1.81188373e+00};
	// c+x+x**2+x**3+x**4
	float pCo[5] = {4.21764218e+00,-9.52098407e-01,5.58119181e-02,-1.62594045e-03,1.74642579e-05};
    // Process the inputs:
    int * inputArrCO = (int *) inputArrCOP;
    int * inputArrA = (int *) inputArrAP;
    int * inputArrADims = (int *) inputArrADimsP;
    float * condArrA = (float *) condArrAP;
    int * inputArrB = (int *) inputArrBP;
    int * inputArrBDims = (int *) inputArrBDimsP;
    float * condArrB = (float *) condArrBP;
    float * outputC = (float *) outputCP;
    
    // set up the variables we will use:
	int xA,posXA,xB,posXB,j,ci;
	int curXAValue,bestDiffxB;
	int minVal,maxVal,lowestI,highestI;
	float posRes,nBlastRes,posD2,nBlastD,posD;
	float compDiff,curDiff,sumDiff2;
	// go through the values in the outData:
	int retval = 0;
	
	for(xA=0;xA<inputArrADims[0];xA++){
		
		posXA = xA*inputArrADims[1];
		
		curXAValue = inputArrA[posXA];
		minVal = curXAValue-threshhold;
		maxVal = curXAValue+threshhold+1;
		minVal = (minVal >= 0) ? minVal : 0;
		maxVal = (maxVal <= xAxis) ? maxVal:xAxis;
		lowestI = inputArrCO[minVal];
		highestI = inputArrCO[maxVal];
		posD2 = 1000000.0;
		bestDiffxB = 0;
		for(xB=lowestI;xB<highestI;xB++){
			posXB = xB*inputArrBDims[1];
			sumDiff2 = 0;
			retval+=1;
			for (j=0;j<3;j++){
				curDiff = inputArrA[posXA+j]-inputArrB[posXB+j];
				curDiff = curDiff*mult[j];
				sumDiff2+=(curDiff*curDiff);
			}
			if (sumDiff2<posD2){
				posD2 = sumDiff2;
				bestDiffxB = xB;
			}
		}
		posD = sqrt(posD2);
		// Only check best Nearest Neighbour!
		if (posD<=posThreshold){
			// c+x+x**2+x**3+x**4
			posRes = pCo[0]+pCo[1]*posD+pCo[2]*pow(posD,2.0)+pCo[3]*pow(posD,3.0)+pCo[4]*pow(posD,4.0);
		}
		else{
			posRes = posScoreIfOverThreshold;
		}
		if (posD<nBlastThreshold){
			nBlastD = 0.0;
			nBlastRes = 0.0;
			for(ci=0;ci<sizeOfComp;ci++){
				compDiff = condArrA[xA*sizeOfComp+ci]*condArrB[bestDiffxB*sizeOfComp+ci];
				nBlastD +=compDiff;
			}
			nBlastD = fabs(nBlastD);
			// c+x+y+x**2+x*y+y*y, x dist squared, y dot product
			nBlastRes = nBCo[0]+nBCo[1]*posD+nBCo[2]*nBlastD+posD*posD*nBCo[3]+posD*nBlastD*nBCo[4]+nBlastD*nBlastD*nBCo[5]+\
			pow(posD,3)*nBCo[6]+pow(posD,2)*nBlastD*nBCo[7]+pow(nBlastD,2)*posD*nBCo[8]+pow(nBlastD,3)*nBCo[9];
		}
		else{
			nBlastRes = nBlastScoreIfOverThreshold;
			nBlastD = -1.0;
		}
		outputC[xA*4] = posD2;
		outputC[xA*4+1]= nBlastD;
		outputC[xA*4+2] = posRes;
		outputC[xA*4+3]= nBlastRes;
	}
	return retval;


}
