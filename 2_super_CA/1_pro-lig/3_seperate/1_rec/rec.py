import os

def Read_PDB(pdb):
  chns=['A','B']
  ress=[8,23,25,27,28,29,30,32,47,48,49,50,81,82,84]
  fr=open("../../2_super/pdb/%s"%(pdb),"r")
  crd=[]
  for line in fr:
    lx=line.split()
    if (lx[0]=='ATOM'):
      if ((lx[4] in chns)&(int(lx[5]) in ress)):
        crd.append(line)
  fr.close()
  return crd

def Write_PDB(crd,pdb):
  fw=open("pdb/%s"%(pdb),"w")
  for line in crd:
    fw.write("%s"%(line))
  fw.close()

def Main():
  global res,pdb,chn
  pdbs=sorted(os.listdir("../../2_super/pdb"))
  for pdb in pdbs:
    crd=Read_PDB(pdb)
    Write_PDB(crd,pdb)

Main()


