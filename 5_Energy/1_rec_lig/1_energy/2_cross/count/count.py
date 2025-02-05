import os
import numpy as np

def read_E():
  fr=open('../Energy')
  rec=[]
  lig=[]
  E=[]
  for line in fr:
    lx=line.split()
    rec.append(lx[0])
    lig.append(lx[1])
    if lx[0]!=lx[1]:
      E.append(float(lx[2]))
  recn=np.array(rec).reshape(79,79)
  En=np.array(E).reshape(79,78)
  return recn,lig,En

def Out_Count(rec,c1,c2):
  fw=open('count','w')
  for i in range (len(c1)):
    fw.write('%4s %6.2f %6.2f\n'%(rec[i][0],c1[i],c2[i]))
  fw.close()

def Main():
  rec,lig,E=read_E()
  count1 = 100*np.sum(E < 0, axis=1)/78
  count2 = 100*np.sum(E > 100, axis=1)/78
  Out_Count(rec,count1,count2)

Main()

