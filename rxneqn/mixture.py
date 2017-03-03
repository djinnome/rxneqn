import copy
import re
from fractions import Fraction
from . import ChemicalFormula
import pandas as pd
class Mixture:
    def __init__( self, mixture ):
        self.mixture = self.parse_mixture( mixture )
    
    def get_species( self ):
        return [m['chemical_formula'] for m in self.mixture]
    def add_to_mixture( self, part ):
        mix = Mixture(str(self))
        for p in part.get_mixture():
            is_new_part = True
            for m in mix.get_mixture():
                if str(m['chemical_formula']) == str(p['chemical_formula']):
                    m['stoichiometry'] += p['stoichiometry']
                    is_new_part = False
            if is_new_part:
                mix.mixture.append( p )
        return mix
    
    def multiply_factor( self, factor ):
        mix = Mixture( str( self ))
        for f in factor.index:
            for m in mix.get_mixture():
                if str(m['chemical_formula']) == str(f):
                    m['stoichiometry'] *= factor[f]
        return mix
    
    def __mul__( self, other ):
        if type( other ) in [int, float, Fraction ]:
            return self.multiply_factor( pd.Series(Fraction( other ), index=self.get_species, dtype='object'))
        if type(other) is pd.Series:
            return self.multiply_factor( other )
        
    def subtract_from_mixture( self, part ):
        mix = Mixture( str( self ))
        for p in part.get_mixture():
            for m in mix.get_mixture():
                if str(m['chemical_formula']) == str( p['chemical_formula']):
                    if m['stoichiometry'] > p['stoichiometry']:
                        m['stoichiometry'] -= p['stoichiometry']
                    else:
                        mix.mixture.remove(  m )
        return mix
                        
                        
    def parse_mixture( self, mixture ):
        if type(mixture) is Mixture:
            mixture = str(mixture)
        pattern = re.compile(r'([0-9./]*\s*)(.*)')
        parts = []
        for part in mixture.split(' + '):
            m = pattern.search(part.strip())
            if m:
                if m.group(1) == '':
                    stoichiometry = Fraction(1)
                else:
                    stoichiometry = Fraction(m.group(1))
                chemical_formula = ChemicalFormula( m.group(2) )
            if stoichiometry != 0:
                parts.append(dict(chemical_formula=chemical_formula, stoichiometry=stoichiometry))
        return parts
    
    def get_mixture( self ):
        return self.mixture
    
    def get_number_of_atoms( self, element ):
        count = 0
        for m in self.mixture:
            count += m['stoichiometry'] * m['chemical_formula'].get_number_of_atoms( element )
        return count
    def get_mass_of_mixture( self ):
        counts = []
        for m in self.mixture:
            counts.append( m['stoichiometry']*m['chemical_formula'].get_molecular_mass() )
        return counts
    
    def get_stoichiometry_of_species(self, species ):
        for m in self.mixture:
            if str(m['chemical_formula']) == str(species):
                return m['stoichiometry']
        return 0
    
    def get_chemical_composition( self ):
        composition = {}
        for m in self.mixture:
            composition[str(m['chemical_formula'])] = m['stoichiometry']*m['chemical_formula'].get_chemical_composition()
        return pd.DataFrame(composition).fillna(0)
    def get_charge_of_mixture( self ):
        charges = 0
        for m in self.mixture:
            charges += m['stoichiometry']*m['chemical_formula'].get_charge()
        return charges

    def to_latex( self ):
        out = []
        for m in sorted(self.mixture, key=lambda x: str(x['chemical_formula'])):
            part = ''
            if m['stoichiometry'] == 1:
                part = str(m['chemical_formula'])
            else:
                part = str(m['stoichiometry']) + '\ ' + str(m['chemical_formula'])
            out.append( part )
        return ' + '.join( out )
    def __repr__( self ):
        out = []
        for m in sorted(self.mixture, key=lambda x: str(x['chemical_formula'])):
            part = ''
            if m['stoichiometry'] == 1:
                part = str(m['chemical_formula'])
            elif m['stoichiometry'] != 0:
                part = str(m['stoichiometry']) + ' ' + str(m['chemical_formula'])
            out.append( part )
        return ' + '.join( out )
    
    def __str__( self ):
        return self.__repr__()
    get_charge = get_charge_of_mixture
    __add__ = add_to_mixture
    __sub__ = subtract_from_mixture
