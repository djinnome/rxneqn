import re
import pandas as pd
import numpy as np
from pint import UnitRegistry
from .utils import LCD
from fractions import Fraction
from . import Mixture, ChemicalFormula
class Reaction:
    u = UnitRegistry()
    u.define('electron-equivalent = [] = e-eq = e-equiv = eeq')
    gas_constant = 8.3144598*u.joule/u.mole/u.kelvin
    def __init__(self, rxn, deltaGEE=None ):
        self.rxn = self.parse_reaction( rxn )
        self.setDeltaGEE( deltaGEE )
    def setDeltaGEE( self, deltaGEE):
        self.deltaGEE = deltaGEE

    def get_deltaGEE_from_potentials( self, mu0 ):
        deltaGEE = 0
        for species in self.get_species():
            try:
                deltaGEE +=  mu0[str(species)]*self.get_stoichiometry_of_species( species )
            except KeyError:
                print("Species {} not in mu0".format( species ))
                deltaGEE += 0
        return deltaGEE
            
    def hasDeltaGEE( self ):
        return self.deltaGEE is not None
    def getDeltaGEE( self ):
        return self.deltaGEE

    # def setDeltaG0( self, deltaG0, temperature=273):
    #     self.deltaGEE = deltaG0 + gas_constant*temperature*u.kelvin*np.ln(10**-7)
    # def getDeltaG0( self, temperature=273):
    #     u = self.u
    #     return self.deltaGEE - self.gas_constant*temperature*u.kelvin*np.ln(10**-7)
        
    def __mul__(self, other):
        if type(other) in [int, float, Fraction]:
            return self.multiply_constant( other )
        elif type(other) is pd.Series:
            return self.multiply_factor( other )
        else:
            raise NotImplemented('unsupported operand type(s) for *: "{}" and "{}"'.format())
    def __rmul__(self, other):
        if type(other) in [int, float, Fraction]:
            return self.multiply_constant( other )
        elif type(other) is pd.Series:
            raise NotImplemented

    def multiply_constant( self, constant ):
        return self.multiply_factor(pd.Series(Fraction(constant),index=self.get_species(),dtype='object'))
    def multiply_factor( self, factor ):
        rxn = Reaction( str(self ) )
        rxn.rxn['reactant'] = rxn.rxn['reactant'].multiply_factor( factor )
        rxn.rxn['product'] = rxn.rxn['product'].multiply_factor( factor )
        return rxn
    
    def add_reaction( self, rxn ):
        new_reactants = Mixture(str(self.rxn['reactant'])).add_to_mixture(rxn.rxn['reactant']).\
            subtract_from_mixture(self.rxn['product']).subtract_from_mixture(rxn.rxn['product'])
        new_products = Mixture( str(self.rxn['product'])).add_to_mixture(rxn.rxn['product']).\
            subtract_from_mixture(self.rxn['reactant']).subtract_from_mixture(rxn.rxn['reactant'])
        if self.hasDeltaGEE() and rxn.hasDeltaGEE():
            newDeltaGEE = self.getDeltaGEE() + rxn.getDeltaGEE()
        else:
            newDeltaGEE = None
        return Reaction(str(new_reactants) + ' ==> ' + str(new_products), newDeltaGEE)
    def subtract_reaction( self, rxn ):
        new_reactants = Mixture(str(self.rxn['reactant'])).add_to_mixture(rxn.rxn['product']).\
            subtract_from_mixture(self.rxn['product']).subtract_from_mixture(rxn.rxn['reactant'])
        new_products = Mixture( str(self.rxn['product'])).add_to_mixture(rxn.rxn['reactant']).\
            subtract_from_mixture(self.rxn['reactant']).subtract_from_mixture( rxn.rxn['product'])
        if self.hasDeltaGEE() and rxn.hasDeltaGEE():
            newDeltaGEE = self.getDeltaGEE() - rxn.getDeltaGEE()
        else:
            newDeltaGEE = None
        return Reaction(str(new_reactants) + ' ==> ' + str(new_products), newDeltaGEE)

    def get_charge( self ):
        return self.rxn['product'].get_charge() - self.rxn['reactant'].get_charge()
    def get_reactants( self ):
        return self.rxn['reactant'].get_species()
    
    def get_products( self ):
        return self.rxn['product'].get_species()
    
    def get_stoichiometry_of_species( self, species ):
        return self.rxn['product'].get_stoichiometry_of_species( species ) - self.rxn['reactant'].get_stoichiometry_of_species(species)
    
    def get_species( self ):
        return self.rxn['reactant'].get_species() + self.rxn['product'].get_species()

    def parse_reaction( self, rxn ):
        if type(rxn) is Reaction:
            rxn = str(Reaction)
        pattern = re.compile(r'==>')
        mixtures = pattern.split(rxn)
        if len(mixtures) == 2:
            reactant = Mixture( mixtures[0].strip())
            product = Mixture( mixtures[1].strip())
            return dict(reactant=reactant, product=product)
    def normalize( self ):
        species = self.get_species()
        lcd = LCD([self.get_stoichiometry_of_species( s ) for s in species])
        return self * lcd
    def deltaGEE_to_latex( self ):
        if self.deltaGEE:
            return r"$\Delta G^{0'} = %0.2f \frac{kJ}{\mathrm{e}^- equiv}$" % self.deltaGEE
        else:
            return ''
    def to_latex( self, deltaGEE=False ):
        if deltaGEE:
            return r'${} \rightarrow {}$ | {}'.format(self.rxn['reactant'].to_latex(), self.rxn['product'].to_latex(), self.deltaGEE_to_latex())
        else:
            return r'${} \rightarrow {}$'.format(self.rxn['reactant'].to_latex(), self.rxn['product'].to_latex()
    def __repr__( self ):
        return str(self.rxn['reactant']) + ' ==> ' + str(self.rxn['product'])
    def __str__( self ):
        return self.__repr__()

    def to_dict( self ):
        rxn = {}
        for species in self.get_species():
            rxn[species] = self.get_stoichiometry_of_species( species )
        return rxn

    def to_series( self ):
        return pd.Series(self.to_dict(), index=self.get_species())

    
    def get_balance( self ):
        return pd.concat( [-self.rxn['reactant'].get_chemical_composition(), 
                           self.rxn['product'].get_chemical_composition()],
                          axis=1).fillna(0).sum(axis=1)

    def get_chemical_composition( self ):
        return pd.concat( [self.rxn['reactant'].get_chemical_composition(), 
                           self.rxn['product'].get_chemical_composition()],
                        axis=1).fillna(0)
    __add__ = add_reaction
    __sub__ = subtract_reaction
