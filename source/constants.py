from rdkit import Chem
# SMARTS representations of the atomic environments for molecule fragmentation.
environs = {
  'L1': '[C;D3]([#0,#6,#7,#8])(=O)',  #original L1
  'L2': '[O;D2]-[#0,#6,#1]',  #original L3
  'L3': '[C;!D1;!$(C=*)]-[#6]',  #original L4
  'L4': '[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]',  #original L5
  'L5': '[C;D2,D3]-[#6]',  #original L7
  'L6': '[C;!D1;!$(C!-*)]',  #original L8
  'L61': '[C;R1;!D1;!$(C!-*)]',
  'L7': '[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]',  #original L9
  'L8': '[N;R;$(N(@C(=O))@[#6,#7,#8,#16])]',  #original L10
  'L9': '[S;D2](-[#0,#6])',  #original L11
  'L10': '[S;D4]([#6,#0])(=O)(=O)',  #original L12
  
  'L11': '[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]',  #original L13
  'L111': '[C;R2;$(C(-;@[C,N,O,S])-;@[N,O,S])]',
  'L112': '[C;R1;$(C(-;@[C,N,O,S;R2])-;@[N,O,S;R2])]',
  
  'L12': '[c;$(c(:[c,n,o,s]):[n,o,s])]',  #original L14
  
  'L13': '[C;$(C(-;@C)-;@C)]',  #original L15
  'L131': '[C;R2;$(C(-;@C)-;@C)]',
  'L132': '[C;R1;$(C(-;@[C;R2])-;@[C;R2])]',
  
  'L14': '[c;$(c(:c):c)]', #original L16
}
reactionDefs = (
  # L1
   [('1', '2', '-'),
    ('1', '4', '-'),
    ('1', '8', '-'),
    ('1', '11', '-'),
    ('1', '12', '-'),
    ('1', '13', '-'),
    ('1', '14', '-')],

  # L2
    [('2', '3', '-'),
    ('2', '11', '-'),
    ('2', '12', '-'),
    ('2', '13', '-'),
    ('2', '14', '-')],

  # L3
    [('3', '4', '-'),
    ('3', '9', '-')],

  # L4
    [('4', '10', '-'),
    ('4', '12', '-'),
    ('4', '14', '-'),
    ('4', '11', '-'),
    ('4', '13', '-')],

  # L5
    [('5', '5', '=')],

  # L6
    [('6', '7', '-'),
    ('6', '8', '-'),
    ('6', '11', '-;!@'),
    ('6', '12', '-'),
    ('6', '13', '-;!@'),
    ('6', '14', '-')],
     
  # L61
    [('61', '111', '-;@'),
    ('61', '131', '-;@')],   

  # L7
    [('7', '11', '-'),  
    ('7', '12', '-'),  
    ('7', '13', '-'),
    ('7', '14', '-')],

  # L8
    [('8', '11', '-'),
    ('8', '12', '-'),
    ('8', '13', '-'),
    ('8', '14', '-')],

  # L9
    [('9', '11', '-'),
    ('9', '12', '-'),
    ('9', '13', '-'),
    ('9', '14', '-')],

  # L11
    [('11', '12', '-'),
    ('11', '13', '-;!@'),
    ('11', '14', '-')],
     
  # L112
    [('112', '132', '-;@')],

  # L12
    [('12', '12', '-'),  
    ('12', '13', '-'),
    ('12', '14', '-')],

  # L13
    [('13', '14', '-')],

  # L14
    [('14', '14', '-')],  
    )

def create_envs_and_bond_matchers():
    environMatchers = {}
    for env, sma in environs.items():
        environMatchers[env] = Chem.MolFromSmarts(sma)

    bondMatchers = []
    for compats in reactionDefs:
        
        tmp = []
        for i1, i2, bType in compats:
            e1 = environs['L%s' % i1]
            e2 = environs['L%s' % i2]
            patt = '[$(%s)]%s[$(%s)]' % (e1, bType, e2)
            patt = Chem.MolFromSmarts(patt)
            tmp.append((i1, i2, bType, patt))
        bondMatchers.append(tmp)
    return environMatchers, bondMatchers


ENVIRONMATCHERS, BONDMATCHERS = create_envs_and_bond_matchers()
