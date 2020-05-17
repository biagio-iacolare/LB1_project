#!/usr/bin/python3

import sys
import math

def get_blast(filename):
    f_list=[]
    d={}
    f=open(filename, 'r')
    for line in f:
        v=line.rstrip().split()
        d[v[0]]=d.get(v[0],[])                              # v[0] = sp|XXXXXX|NAME_
        d[v[0]].append([float(v[1]),int(v[2])])             # d[v[0]] = [best-E-value, 0/1]
    for k in d.keys():
        d[k].sort()
        f_list.append(d[k][0])
    return f_list

def get_cm(data,th):        
    cm=[[0.0,0.0],[0.0,0.0]]                    # CM= Confusion Matrix = [[TP,FP],[FN,TN]]
                                                # 1 and 0 are the two classes of positives and negatives
    for i in data:
        if i[0]<th and i[1]==1:                 # E-value<th & class=1 (Kunitz)      = is a TP
            cm[0][0]+=1
        if i[0]>=th and i[1]==1:                # E-value>=th & class=1 (Kunitz)     = is a FN
            cm[1][0]+=1
        if i[0]<th and i[1]==0:                 # E-value<th & class=0 (NOT Kunitz)  = is a FP
            cm[0][1]+=1
        if i[0]>=th and i[1]==0:                # E-value>=th & class=0 (NOT Kunitz) = is a TN
            cm[1][1]+=1
    return cm

def get_acc(cm):
    return float(cm[0][0]+cm[1][1])/(sum(cm[0])+sum(cm[1]))                     # Acc = (TP+TN) / (totP + totN)

def mcc(cm):
  d=(cm[0][0]+cm[1][0])*(cm[0][0]+cm[0][1])*(cm[1][1]+cm[1][0])*(cm[1][1]+cm[0][1])     # d = denominator = (TP+FN)*(TP+FP)*(TN+FN)*(TN+FP)
  return (cm[0][0]*cm[1][1]-cm[0][1]*cm[1][0])/math.sqrt(d)                             # MCC = (TP*TN - FP*FN) / sqrt(d)


if __name__=='__main__':
    filename=sys.argv[1]
    data=get_blast(filename)                                                # data is a list of 2-elements-lists, [ E-value , class ] for each protein
    for i in range(20):
        th=10**-i                                                           # test all thresholds from 1 to 1e-20
        cm=get_cm(data,th)
        print('TH:',th,'\t','Acc:',get_acc(cm),'\t','MCC:',mcc(cm),'\t',cm) 
