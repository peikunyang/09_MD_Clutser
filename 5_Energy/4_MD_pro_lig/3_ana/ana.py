import os
import numpy as np

rec=[0.4,0.6,0.8,1.0,1.2,1.4]
lig=[0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]

def Read_cluster(fn):
  E=[]
  fr=open(fn,'r')
  for line in fr:
    lx=line.split()
    E.append(float(lx[3]))
  fr.close()
  return E

def Write_Data(num,E):
  fw=open('Energy/Energy_%3.1f'%(rec[num]),'w')
  for j in range (len(lig)):
    fw.write('   lig_%3.1f'%(lig[j]))
  fw.write('\n')
  for i in range (E.shape[1]):
    for j in range (E.shape[0]):
      fw.write(' %9.2f'%(E[j][i]))
    fw.write('\n')
  fw.close()

def main():
  for i in range (len(rec)):
    En=[]
    for j in range (len(lig)):
      E=Read_cluster('../2_clu_min/rec_%3.1f/lig_%3.1f/Energy'%(rec[i],lig[j]))
      En.append(E)
    En2=np.array(En)
    Write_Data(i,En2)

main()

