import os
import numpy as np

def Read_1_pro_lig():
  pdb=[]
  E=[]
  fr=open('../../1_rec_lig/1_energy/1_self/Energy','r')
  for line in fr:
    lx=line.split()
    pdb.append(lx[0])
    E.append(float(lx[1]))
  fr.close()
  return pdb,E

def Read_1_all_trj():
  pdb=[]
  E=[]
  fr=open('../1_all_trj/energy','r')
  for line in fr:
    lx=line.split()
    pdb.append(lx[1])
    E.append(float(lx[3]))
  fr.close()
  return pdb,E

def Read_2_cluster(cluster):
  pdb=[]
  E=[]
  fr=open('../2_cluster/%s/Energy'%(cluster),'r')
  for line in fr:
    lx=line.split()
    pdb.append(lx[0][:4])
    E.append(float(lx[2]))
  fr.close()
  return pdb,E

def Chk(pdb):
  for i in range (pdb.shape[1]):
    for j in range (1,pdb.shape[0]):
      if pdb[j][i]!=pdb[0][i]:
        print(pdb[j][i])

def Write_Data(pdb,E):
  fw=open('Energy','w')
  fw.write('          Ori      all      0.4      0.5      0.6      0.7      0.8      0.9      1.0      1.1      1.2      1.3      1.4      1.5      1.6\n')
  for i in range (pdb.shape[0]):
    fw.write('%4s'%(pdb[i][:4]))
    for j in range (E.shape[0]):
      fw.write(' %8.2f'%(E[j][i]))
    fw.write('\n')
  fw.close()

def main():
  cluster=['div_0.4','div_0.5','div_0.6','div_0.7','div_0.8','div_0.9','div_1.0','div_1.1','div_1.2','div_1.3','div_1.4','div_1.5','div_1.6']
  pdbn=[]
  En=[]

  pdb,E=Read_1_pro_lig()
  pdbn.append(pdb)
  En.append(E)

  pdb,E=Read_1_all_trj()
  pdbn.append(pdb)
  En.append(E)

  for i in range (len(cluster)):
    pdb,E=Read_2_cluster('%s'%(cluster[i]))
    pdbn.append(pdb)
    En.append(E)

  pdbn2=np.array(pdbn)
  En2=np.array(En)
  Chk(pdbn2)
  Write_Data(pdbn2[0],En2)

main()

