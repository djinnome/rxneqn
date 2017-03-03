from rxneqn import HalfReactionBalancer, Reaction
import re
import pandas as pd

alanine = 'CH3CHNH2COOH'
carbon_dioxide = 'CO2'
nitrogen_source = 'NH3'
half_rxn = HalfReactionBalancer() #carbon_dioxide, alanine, 'C', nitrogen_source)
truth_eqn = """1/24 C6H12O6 + 1/4 H2O ==> 1/4 CO2 + E- + H+
E- + 6/5 H+ + 1/5 NO3- ==> 3/5 H2O + 1/10 N2
1/24 C6H12O6 + 1/5 H+ + 1/5 NO3- ==> 1/4 CO2 + 7/20 H2O + 1/10 N2
1/4 CO2 + E- + H+ + 1/12 NH3 ==> 1/12 CH3CHNH2COOH + 1/3 H2O
1/5 CO2 + E- + H+ + 1/20 HCO3- + 1/20 NH4+ ==> 1/20 C5H7O2N + 9/20 H2O
1/8 CO2 + E- + H+ ==> 1/8 CH4 + 1/4 H2O
1/6 CO2 + E- + H+ + 1/12 HCO3- + 1/12 NH4+ ==> 1/12 CH3CHNH2COOH + 5/12 H2O
1/12 CH3CHNH2COOH + 1/6 H2O ==> 1/8 CH4 + 1/24 CO2 + 1/12 HCO3- + 1/12 NH4+
1/5 CO2 + E- + H+ + 1/10 HCO3- ==> 1/10 CH3COCOO- + 2/5 H2O
1/4 CO2 + E- + H+ + 1/12 NH3 ==> 1/12 CH3CHNH2COOH + 1/3 H2O
1/6 CO2 + E- + H+ + 1/12 HCO3- + 1/12 NH4+ ==> 1/12 CH3CHNH2COOH + 5/12 H2O
1/12 CO2 + E- + H+ + 1/6 HCO3- + 1/12 NH4+ ==> 1/12 CH3CHNH2COO- + 1/2 H2O
1/5 CO2 + E- + H+ + 1/20 HCO3- + 1/20 NH4+ ==> 1/20 C5H7O2N + 9/20 H2O
4/19 CO2 + E- + 20/19 H+ + 1/19 HCO3- + 1/19 NH4+ ==> 1/19 C5H7O2N+ + 9/19 H2O
""".split('\n')
eqn = {}
eqn[13] = half_rxn.balance_half_reaction('C6H12O6', 'CO2')
eqn[14] = half_rxn.balance_half_reaction('NO3-', 'N2',nitrogen_source='NO3-')
eqn[15] = eqn[13] +  eqn[14] 
eqn[16] = half_rxn.balance_half_reaction(carbon_dioxide, alanine, nitrogen_source)
eqn[17] = half_rxn.balance_half_reaction('CO2', 'C5H7O2N','NH4+')
table3 = {'O-12': '1/8 CO2 + H+ + E- ==> 1/8 CH4 + 1/4 H2O',
         'O-2': '1/6 CO2 + 1/12 NH4+ + 1/12 HCO3- + H+ + E- ==> 1/12 CH3CHNH2COOH + 5/12 H2O',
          'O-3': str(half_rxn.balance_half_reaction('CO2','C6H5COO-')),
          'O-20': '1/5 CO2 + E- + H+ + 1/20 HCO3- + 1/20 NH4+ ==> 1/20 C5H7O2N + 9/20 H2O'
         }

eqn[18] = Reaction(table3['O-12'])
eqn[19] = Reaction(table3['O-2'])
eqn[20] = eqn[18] - eqn[19]
eqn[21] = half_rxn.balance_half_reaction('CO2','CH3COCOO-')
eqn[22] = half_rxn.balance_half_reaction('CO2','CH3CHNH2COOH','NH3')
eqn[23] = half_rxn.balance_half_reaction('CO2','CH3CHNH2COOH','NH4+')
eqn[24] = half_rxn.balance_half_reaction('CO2','CH3CHNH2COO-','NH4+')
eqn[25] = half_rxn.custom_half_reaction( C=5,H=7, O=2, N=1)
eqn[26] = half_rxn.custom_half_reaction( C=5,H=7, O=2, N=1, charge=1)
k = 0

molecules = dict(glucose='C6H1206',
                NAD='C21H26N7O14P2',
                ADP='C10H12N5O10P2-3',
                Pi='HO4P-2',
                proton='H+',
                pyruvate='C3H3O3',
                NADH='C21H27N7O14P2',
                ATP='C10H12N5O13P3-4',
                water='H2O')
atp_hydrolysis = Reaction('{ATP} + {water} ==> {ADP} + {Pi} + {proton}'.format(**molecules))
balance = pd.Series(dict(
    C      =   0.0,
    Charge =   0.0,
    H      =   0.0,
    N      =   0.0,
    O      =   0.0,
    P      =   0.0))

for i in eqn:
    assert str(eqn[i]) == truth_eqn[k], "Equation[{}]: {} should be {}".format(i,str(eqn[i]), truth_eqn[k])
    k += 1
assert(atp_hydrolysis.get_balance() == balance, "{} should be {}".format(atp_hydrolysis.get_balance(), balance))
print("{} tests out of {} passed".format(k+1, len(eqn)+1))
