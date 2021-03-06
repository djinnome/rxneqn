__version__= '0.1'
import re
import sys
import pandas as pd
from fractions import Fraction, gcd

import periodic
from periodic import element
from .chemical_formula import ChemicalFormula
from .mixture import Mixture
from .reaction import Reaction
from .utils import LCD, GCD, LCM
from .half_rxn_balancer import HalfReactionBalancer


molecular_formula = dict(glucose='C6H12O6',
                         ammonia='NH3',
                         ammonium='NH4+',
                         bicarbonate='HCO3-',
                         oxygen='O2',
                         sucrose='C12H22O11',
                         NAD='C21H26N7O14P2',
                         ADP='C10H12N5O10P2-3',
                         Pi='H2O4P-',
                         carbon_dioxide = 'CO2',
                         proton='H+',
                         NADH='C21H27N7O14P2',
                         ATP='C10H12N5O13P3-4',
                         water='H2O',
                         acetate='CH3COO-',
                         alanine='CH3CHNH2COO-',
                         benzoate='C6H5COO-',
                         citrate='C6H5O7-3',
                         ethanol='CH3CH2OH',
                         formate='HCOO-',
                         glutamate='COOHCH2CH2CHNH2COO-',
                         glycerol='CH2OHCHOHCH2OH',
                         glycine='CH2NH2COOH',
                         lactate='CH3CHOHCOO-',
                         methane='CH4',
                         methanol='CH3OH',
                         palmitate='C16H31O2-',
                         propionate='CH3CH2COO-',
                         pyruvate='CH3COCOO-',
                         succinate='C4H4O4-2',
                         nitrate='NO3-',
                         sulfate='SO4-2',
                         ironIII='Fe+3',
                         ironII='Fe+2',
                         hydrogen='H2',
                         dinitrogen='N2')

mu0m = {'ADP': -1422.5, # Equilibrator mu0'm (ph7, 1 mM concentrations)
 'ATP': -2295.1,
 'Pi': -1052.8,
 'ammonia': -26.5,
 'ammonium': -79.0,
 'carbon_dioxide': -403.1,
 'glucose': -446.8,
 'proton': -17.1,
 'water': -157.6}


mu0ph7 = dict(Pi=-1052.8, # Equilibrator mu0' (ph7, 1 M concentrations)
              ATP=-2295.8,
              ADP=-1423.6,
              glucose=-429.7,
              ammonia=89.6,
              ammonium=79.0,
              water=-157.6,
              carbon_dioxide=-386.0,
              bicarbonate=-546.8,
              proton=0,
              hydrogen=98.7
)

           

deltaGEE = dict(
                glucose=41.35,
                acetate=27.40,
                alanine=31.37,
                benzoate=27.34,
                citrate=33.08,
                ethanol=31.18,
                formate=39.19,
                ammonia=-26.57,
                ammonium=-79.37,
                glutamate=30.93,
                glycerol=38.88,
                glycine=39.80,
                lactate=32.29,
                methane=23.53,
                methanol=36.84,
                palmitate=27.26,

                propionate=27.63,
                pyruvate=35.09,
                succinate=29.09,
                oxygen=-78.72,
                nitrate=-72.20,
                sulfate=20.85,
                iron=-74.27,
                pyruvate_to_cell_using_ammonium=18.8
               )
# def make_half_reactions(molecules, deltaGEE):
#     half_rxn = HalfReactionBalancer()
#     rxns = {}
#     for molecule in molecules:
#         rxns[molecule] = half_rxn.balance_half_reaction(oxidized_form='CO2',reduced_form=molecules[molecule])
#         if molecule in deltaGEE:
#             rxns[molecule].setDeltaGEE(deltaGEE[molecule])
#     return rxns

# half_reactions = make_half_reactions( molecular_formula, deltaGEE )
