
# MacFrag is an efficient molecule fragmentation method, which is capable of segmenting large-scale molecules in a rapid speed and 
# generating diverse fragments that are more compliant with the Rule of Three.

# MacFrag is developed and maintained by Prof. HongLin Li's Group, School of Pharmacy, East China University of Science & Technology, Shanghai 200237, China. 
# http://lilab-ecust.cn/

import argparse
from .source.fragmentation import write_file

#RDLogger.DisableLog('rdApp.*')


def parse_args():
    parser = argparse.ArgumentParser(description='Molecular Fragmentation')
    parser.add_argument('-input_file', '-i', required=True, type=str,                       
                        help='.smi or .sdf file of molecules to be fragmented')    
    
    parser.add_argument('-output_path', '-o', required=True, type=str,
                        help='path of the output fragments file')
    
    parser.add_argument('-maxBlocks', default=1, type=int,
                        help='the maximum number of building blocks that the fragments contain')
    
    parser.add_argument('-maxSR', default=8, type=int,
                        help='only cyclic bonds in smallest SSSR ring of size larger than this value will be cleaved')

    parser.add_argument('-asMols',  default=True, type=bool,
                        help='True of False; if True, MacFrag will reture fragments as molecules and the fragments.sdf file will be output; if False, MacFrag will reture fragments.smi file with fragments representd as SMILES strings')
    
    parser.add_argument('-minFragAtoms',  default=3, type=int,
                        help='the minimum number of atoms that the fragments contain')

    return parser.parse_args()


if __name__ == '__main__':    
    args = parse_args()
    write_file(args.input_file, args.output_path,
               args.maxBlocks, args.maxSR, args.asMols, args.minFragAtoms)
