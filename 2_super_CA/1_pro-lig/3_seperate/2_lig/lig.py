import os

def Read_Lig():
  fr=open("ligand_info","r")
  pdbs=[]
  ligs=[]
  for line in fr:
    lx=line.split()
    pdbs.append(lx[0])
    num=int((len(lx)-1)/3)
    lig=[]
    for i in range (num):
      lig.append((lx[3*i+1],lx[3*i+2]))
    ligs.append(lig)
  fr.close()
  return pdbs,ligs

def Read_PDB(pdb):
  fr=open("../../2_super/pdb/%s.pdb"%(pdb),"r")
  crd=[]
  for line in fr:
    if ((line[:4]=='ATOM')|(line[:6]=='HETATM')):
      if ((line[21]!='A')&(line[21]!='B')):
        crd.append(line)
  fr.close()
  return crd

def Write_PDB(crd,pdb,lig):
  for i in range (len(lig)):
    chn=lig[i][0]
    res=lig[i][1]
    fw=open("pdb/%s.pdb"%(pdb[:4]),"w")
    for line in crd:
      if ((line[21]==chn)&(line[17:20]==res)):
        fw.write("%s"%(line))
    fw.close()

def Main():
  global res,pdb,chn
  pdb,lig=Read_Lig()
  for i in range (len(pdb)): # len(pdb)
    crd=Read_PDB(pdb[i])
    Write_PDB(crd,pdb[i],lig[i])

Main()

