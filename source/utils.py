from rdkit import Chem
from itertools import combinations
import copy

from .constants import *


def SSSRsize_filter(bond,maxSR=8):
    
    judge=True
    for i in range(3,maxSR+1):
        if bond.IsInRingSize(i) :
            judge=False
            break           
    return judge
       

def searchBonds(mol,maxSR=8):
    
    bondsDone = set()
    
    envMatches = {}   
    for env, patt in ENVIRONMATCHERS.items():
        envMatches[env] = mol.HasSubstructMatch(patt)
        
    for compats in BONDMATCHERS:
                
        for i1, i2, bType, patt in compats:
            if not envMatches['L' + i1] or not envMatches['L' + i2]:
                continue
            
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                if match not in bondsDone and (match[1], match[0]) not in bondsDone:
                    bond=mol.GetBondBetweenAtoms(match[0],match[1])
                    
                    if not bond.IsInRing():
                        bondsDone.add(match) 
                        yield (((match[0], match[1]), (i1, i2)))
                    elif bond.IsInRing() and  SSSRsize_filter(bond,maxSR=maxSR):
                        bondsDone.add(match)
                        yield (((match[0], match[1]), (i1, i2)))
                    
def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def mol_remove_atom_mapnumber(mol):
    for atom in mol.GetAtoms():
        atom.ClearProp('molAtomMapNumber')
    return mol

def Get_block_index(blocks):   
    
    block_index={}
    i=0
    for bs in blocks:
        tmp=[a.GetAtomMapNum() for a in bs.GetAtoms() if a.GetSymbol()!='*']
#        tmp=list(set(tmp))
        block_index[tuple(tmp)]=i
        i+=1      
    return block_index

def extrac_submol(mol,atomList,link):
    aList_mol=list(range(mol.GetNumAtoms()))
    aList_link=list(set([a[0] for a in link]))
    aList_submol=list(set(atomList+aList_link))
    aList_remove=[a for a in aList_mol if a not in aList_submol]
    eMol=Chem.RWMol(mol)
    
    aList_bbond=[a for a in aList_link if a not in atomList]
    for b in combinations(aList_bbond,2):
        eMol.RemoveBond(b[0],b[1])
    
    aList_remove.sort(reverse=True)
    for a in aList_remove:
        eMol.RemoveAtom(a) 

            
    for ba,btype in link:
        if ba in atomList:
            continue
        tmpatom=[a for a in eMol.GetAtoms() if a.GetAtomMapNum()==ba][0]
        tmpatom.SetIsAromatic(False)
        tmpatom.SetAtomicNum(0)
        tmpatom.SetIsotope(int(btype))
        tmpatom.SetNoImplicit(True)
            
    frag=eMol.GetMol()
   
    for a in frag.GetAtoms():
        a.ClearProp('molAtomMapNumber')   
    return frag
 
    
def simple_iter(graph, k):
    cis=[]
#    numb_at_most_k = 0
    # colors for vertices (igraph is to slow)
    colors = graph.vcount()*[-1]
    # consider only graphs with at least k vertices
    for i in range(graph.vcount() -1, -1, -1):
        # subgraph_set
        subgraph_set = [i]
        # print first induced subgraph conisting only of vertex i
#        subgraph_file.write(print_names(graph, subgraph_set))
        cis.append(copy.deepcopy(subgraph_set))
#        numb_at_most_k += 1
        # extension set
        extension = []
        # list of lists of all exclusive neighbors
        ex_neighs = [[]]
        # color closed neighborhood of start vertex
        colors[i] = 0
        for vertex in graph.neighbors(i):
            if vertex >= i:
                break
            colors[vertex] = 1
            extension.append(vertex)
            ex_neighs[0].append(vertex)
        # lists for the pointers of the extension set
        pointers = [extension.__len__()-1]
        # we save the index in whose closed neighborhood the actual pointers are, to perform the jumps
        poi_col = [1]
        # subgraph size
        sub_size = 1
        # enumerate all subgraphs containing vertex i
        while subgraph_set != []:
            # for each vertex in the extension set, create a new branch
            # if the actual pointer points to null (-1), then the corresponding extension set is empty
            while pointers[-1] > -1:
                last = pointers[-1]
                vertex = extension[last]
                ver_col = colors[vertex]
                # move pointer one to the left
                pointers[-1] -= 1
                act_vertex = pointers[-1]
                # check if vertex is an exclusive neighbor of the vertex last added to the subgraph set
                if ver_col == sub_size:
                    # if yes delete this vertex from the extension set
                    extension.pop()
                # check if pointer jumps, in other words check if the new next vertex in the extension set has a different color than the current vertex
                if act_vertex > -1:
                    if colors[extension[act_vertex]] < poi_col[-1]:
                        # jump with pointer (we stay at the same point if the color of the actual vertex is one less)
                        pointers[-1] = pointers[poi_col[-1]-2]
                        # update color
                        poi_col[-1] = colors[extension[pointers[-1]]]
                # create next child
                subgraph_set.append(vertex)
#                numb_at_most_k += 1
#                subgraph_file.write(print_names(graph, subgraph_set))
                cis.append(copy.deepcopy(subgraph_set))
                sub_size += 1
                # check if adding this vertex leads to a solution
                if sub_size == k:
                    # check time limit
#                    time_now = clock()
#                    if time_now > time_max:
 #                       return numb_at_most_k
                    subgraph_set.pop()
                    sub_size -= 1
                # otherwise create new enumeration tree node
                else:
                    # find the exclusive neighbors of vertex
                    ex_neighs.append([])
                    found = False
                    for neig in graph.neighbors(vertex):
                        if neig >= i:
                            break
                        if colors[neig] == -1:
                            colors[neig] = sub_size
                            ex_neighs[-1].append(neig)
                            extension.append(neig)
                            found = True
                    # if there is an exclusive neighbor, the new pointer points to him, otherwise to the old pointer
                    if found == True:
                        pointers.append(extension.__len__()-1)
                        poi_col.append(sub_size)
                    else:
                        pointers.append(pointers[-1])
                        poi_col.append(poi_col[-1])
            # now the last pointer points to null, so restore the parent
            pointers.pop()
            poi_col.pop()
            subgraph_set.pop()
            for vertex in ex_neighs[-1]:
                colors[vertex] = -1
            ex_neighs.pop()
            sub_size -= 1
        # remove color from the root
        colors[i] = -1
    return cis    