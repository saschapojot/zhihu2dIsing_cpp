import xml.etree.ElementTree as ET
import numpy as np
# from datetime import datetime
import statsmodels.api as sm
from scipy import stats
import glob
import sys
import re
from copy import deepcopy
import warnings
# import matplotlib.pyplot as plt

#This script checks if a vector reaches equilibrium

sigWrongArg="wrong number of arguments"
sigEq="equilibrium"
sigContinue="continue"
sigStop="stop"

if (len(sys.argv)!=2):
    print(sigWrongArg)
    exit()

xmlFilesPath=str(sys.argv[1])

#fetch files in the directory
inXMLFileNames=[]
startVals=[]
for file in glob.glob(xmlFilesPath+"/*"):
    inXMLFileNames.append(file)
    matchStart=re.search(r"loopStart(-?\d+(\.\d+)?)loopEnd",file)
    if matchStart:
        startVals.append(matchStart.group(1))


def str2int(valList):
    ret = [int(strTmp) for strTmp in valList]
    return ret


startVals = str2int(startVals)

start_inds = np.argsort(startVals)

#sort files by starting value
inXMLFileNames=[inXMLFileNames[ind] for ind in start_inds]


#ensure the file number is a multiple of 3
if len(inXMLFileNames)%3==0:
    xmlFileToBeParsed=deepcopy(inXMLFileNames)
elif len(inXMLFileNames)%3==1:
    xmlFileToBeParsed=deepcopy(inXMLFileNames[1:])
else:
    xmlFileToBeParsed=deepcopy(inXMLFileNames[2:])

xmlFileToBeParsed=xmlFileToBeParsed[int(len(xmlFileToBeParsed)/3*2):]

# print("xml file number: "+str(len(xmlFileToBeParsed)))
def parse1File(fileName):
    """

    :param fileName: xml file name
    :return: the values in the vector
    """

    tree=ET.parse(fileName)
    root = tree.getroot()
    vec=root.find("vec")
    vec_items=vec.findall('item')
    vecValsAll=[float(item.text) for item in vec_items]
    # vecValsAll=np.array(vecValsAll)
    return vecValsAll

#combine all vectors
vecValsCombined=parse1File(xmlFileToBeParsed[0])

for file in xmlFileToBeParsed[1:]:
    vecValsCombined+=parse1File(file)



vecValsCombined=np.array(vecValsCombined)
# print(len(vecValsCombined))

#ferromagnetic case: all values equal

# meanE=np.mean(vecValsCombined)
# print("mean="+str(meanE))
diff=np.linalg.norm(vecValsCombined-1,ord=1)/len(vecValsCombined)
# print("diff="+str(diff))

if diff<1e-7:
    print(sigStop+" ferro")
    exit()




#check if the whole vector has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        vecAutc=sm.tsa.acf(vecValsCombined)
    except Warning as w:
        print(sigStop+" ferro")
        exit()



halfLength=int(len(vecValsCombined)/2)

part0=vecValsCombined[:halfLength]
part1=vecValsCombined[halfLength:]

ferro0=False
ferro1=False
#check if the part0 has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        autc0=sm.tsa.acf(part0)
    except Warning as w:
        ferro0=True
#check if the part1 has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        autc1=sm.tsa.acf(part1)
    except Warning as w:
        ferro1=True

if ferro0 and ferro1:
    print(sigStop+" ferro")
    exit()
elif ferro0==True and ferro1==False:
    # print("f0 True f1 False")
    print(sigContinue)
    exit()
elif ferro0==False and ferro1==True:
    # print("f0 False f1 True")
    print(sigContinue)
    exit()

def Jackknife(vec):
    """

    :param vec:
    :return: the mean and half length  of 0.95 confidence interval computed by Jackkknife
    """

    psMean=np.mean(vec)

    psVar=np.var(vec,ddof=1)

    n=len(vec)

    hfLen=1.96*np.sqrt(psVar/n)
    return psMean,hfLen


#computation of auto-correlation

# M=100
# lags=30000
NLags=int(len(vecValsCombined)*2/3)
acfOfVec=sm.tsa.acf(vecValsCombined,nlags=NLags)
# print("min correlation is ",np.min(acfOfVec))
# print("total elem number = "+str(len(vecValsCombined)))
# plt.figure()
# plt.plot(acfOfVec,color="black")
# plt.savefig("ECorr.png")
# plt.close()
eps=(1e-2)*5
pThreshHold=0.05
lagVal=0
if np.min(np.abs(acfOfVec))>eps:
    print("high correlation")
    print(np.min(np.abs(acfOfVec)))

    print(sigContinue)
    exit()
else:
    lagVal=np.where(np.abs(acfOfVec)<=eps)[0][0]
    selectedFromPart0=part0[::lagVal]
    # print(len(selectedFromPart0))
    selectedFromPart1=part1[::lagVal]
    # plt.subplot(1,2,1)
    # plt.hist(selectedFromPart0,bins=100)
    # plt.subplot(1,2,2)
    # plt.hist(selectedFromPart1,bins=100)
    # print(np.mean(selectedFromPart0), np.sqrt(np.var(selectedFromPart0))/np.sqrt(len(selectedFromPart0/60)))
    # print(np.mean(selectedFromPart1), np.sqrt(np.var(selectedFromPart1))/np.sqrt(len(selectedFromPart1/60)))

    # plt.savefig("hist.png")
    # D,p=stats.ks_2samp(part0,part1)
    mean0,hf0=Jackknife(part0)
    mean1,hf1=Jackknife(part1)
    print("mean0="+str(mean0)+", mean1="+str(mean1))
    print("hf0="+str(hf0)+", hf1="+str(hf1))
    if np.abs(mean0-mean1)<=hf0 or np.abs(mean0-mean1)<=hf1:
        print(sigEq+" "+str(lagVal))
        exit()
    else:
        print(sigContinue)
        exit()



