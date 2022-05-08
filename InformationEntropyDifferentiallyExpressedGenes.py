from math import log

def calcShannonEnt(dataSet):
    shannonEnt = 0.0
    total = sum(dataSet)
    for k in dataSet:
        prob = k/total
        shannonEnt -= prob * log(prob,2)
    return shannonEnt

dataSet = [8,8,8,8,8,8,8,8]
print(log(8,2) - calcShannonEnt(dataSet))
dataSet = [1,2,3,4,5,6,7,8]
print(log(8,2) - calcShannonEnt(dataSet))
dataSet = [1,1,1,1,1,1,1,8]
print(log(8,2) - calcShannonEnt(dataSet))