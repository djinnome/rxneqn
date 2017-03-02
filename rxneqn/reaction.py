import re
import pandas as pd
import numpy as np
from pint import UnitRegistry
from .utils import LCD
from fractions import Fraction
from . import Mixture, ChemicalFormula
class Reaction:
    def __init__(self, rxn ):
        self.rxn = self.parse_reaction( rxn )

    def setGibbsFreeEnergy( self, deltaG0pH7 ):
        self.unit_registry = UnitRegistry()
        self.deltaG0pH7 = deltaG0pH7
        
    def getStandardFreeEnergy( self, temperature=273):
        u = self.unit_registry
        gas_constant = 8.3144598*u.joule/u.mole/u.kelvin
        return self.deltaG0pH7 - gas_constant*temperature*u.kelvin*np.ln(10**-7)
        
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
        return Reaction(str(new_reactants) + ' ==> ' + str(new_products))
    def subtract_reaction( self, rxn ):
        new_reactants = Mixture(str(self.rxn['reactant'])).add_to_mixture(rxn.rxn['product']).\
            subtract_from_mixture(self.rxn['product']).subtract_from_mixture(rxn.rxn['reactant'])
        new_products = Mixture( str(self.rxn['product'])).add_to_mixture(rxn.rxn['reactant']).\
            subtract_from_mixture(self.rxn['reactant']).subtract_from_mixture( rxn.rxn['product'])
        return Reaction(str(new_reactants) + ' ==> ' + str(new_products))

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
    
    def __repr__( self ):
        return str(self.rxn['reactant']) + ' ==> ' + str(self.rxn['product'])
    def __str__( self ):
        return self.__repr__()

    def get_chemical_composition( self ):
        return pd.concat( [self.rxn['reactant'].get_chemical_composition(), 
                           self.rxn['product'].get_chemical_composition()],
                        axis=1).fillna(0)
    __add__ = add_reaction
    __sub__ = subtract_reaction
