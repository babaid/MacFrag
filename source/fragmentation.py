from click import FileError
from rdkit import Chem
from rdkit.Chem.BRICS import BreakBRICSBonds
from igraph import Graph
import copy

from .utils import *
from .constants import *


def MacFrag(mol, maxBlocks=6, maxSR=8, asMols=True, minFragAtoms=1):
    
    fragPool={}
    
    mol=mol_with_atom_index(mol)    
   
    bonds=list(searchBonds(mol,maxSR=maxSR))
    fragmol=BreakBRICSBonds(mol,bonds=bonds)
    blocks=Chem.GetMolFrags(fragmol,asMols=True)
    for block in blocks:
        tmp=copy.deepcopy(block)
        tmp_smiles=Chem.MolToSmiles(tmp)  
        nAtoms=tmp.GetNumAtoms(onlyExplicit=True)
        if nAtoms - tmp_smiles.count('*') >= minFragAtoms:
            for a in tmp.GetAtoms():
                a.ClearProp('molAtomMapNumber') 
            fragPool[tmp]=Chem.MolToSmiles(tmp)

        
    if maxBlocks>len(blocks):
        maxBlocks=len(blocks)
    if maxBlocks==1:
        if asMols:
            return fragPool.keys()
        else:
            return list(set(fragPool.values()))
       
    block_index=Get_block_index(blocks)
    index_block={block_index[key]:key for key in block_index.keys()}
                

    bond_block={}
    point=False
    block_link={block:set() for block in block_index.keys()}
    for b in bonds:
        ba1,ba2=b[0]
        btype1,btype2=b[1]
        for block in block_index.keys():
            if ba1 in block:
                bond_block[ba1]=block
                block_link[block].add((ba2,btype2))
                point=True
            if ba2 in block:
                bond_block[ba2]=block
                block_link[block].add((ba1,btype1))
                point=True
            if point==True:continue
    
    n=len(index_block.keys())    
    edges=[]
    for b in bonds:
        ba1,ba2=b[0]
        edges.append((block_index[bond_block[ba1]],block_index[bond_block[ba2]]))
    
    graph=Graph(n=n,edges=edges,directed=False)
#    all_cis=rssp(graph,k=maxBlocks)
    all_cis=simple_iter(graph,k=maxBlocks)
    for i in range(n):
        all_cis.remove([i])
    
    sub_link={}
    for cis in all_cis:
        tmp=[]
        tmp2=set()
        for ni in cis:
            tmp.extend(list(index_block[ni]))
            tmp2=tmp2.union(block_link[index_block[ni]])
        sub_link[tuple(tmp)]=tmp2  

  
    for fa in sub_link:
        frag=extrac_submol(mol,atomList=list(fa),link=sub_link[fa])
        frag_smiles=Chem.MolToSmiles(frag)    
        nAtoms=frag.GetNumAtoms(onlyExplicit=True)
        if nAtoms - frag_smiles.count('*') >= minFragAtoms:
            fragPool[frag]=frag_smiles
                
        
    if asMols:
        return list(fragPool.keys())
    else:
        return list(set(fragPool.values()))

def write_file(input_file, dir, maxBlocks, maxSR, asMols, minFragAtoms):

    if input_file.endswith(".smi"):
        smiles=[line.strip().split()[0] for line in open(input_file).readlines()]
        mols=[Chem.MolFromSmiles(s) for s in smiles]
    elif input_file.endswith("mol2"):
        mols = [Chem.MolFromMol2File(input_file, sanitize=False)]
    elif input_file.endswith("sdf"):
        mols = Chem.SDMolSupplier(input_file)
    else:
        raise FileError(f"The extension of the provided file ({input_file}) is wrong. Possible extensions: SDF or MOL2 OR SMI")
    
    for i, mol in enumerate(mols):
        if mol is None:
            print(f"Molecule {i} is not loaded properly.")

    frags = None
    if asMols == True:
        fw = Chem.SDWriter(dir)
        for mol in mols:
            if mol is not None:
                frags = MacFrag(mol, maxBlocks=maxBlocks, maxSR=maxSR, asMols=True, minFragAtoms=minFragAtoms)
            if frags is not None:
                for f in frags:
                    fw.write(f)
    else:
        for mol in mols:
            if mol is not None:
                frags = MacFrag(mol, maxBlocks=maxBlocks, maxSR=maxSR, asMols=False, minFragAtoms=minFragAtoms)
        if frags is not None:
            with open(dir, "w") as file:
                for f in frags:
                    file.write(f+'\n')
        