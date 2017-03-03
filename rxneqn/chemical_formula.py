import re
import pandas as pd
from fractions import Fraction
import periodic
from periodic import element

class ChemicalFormula:
    chemical_formula_RE = r'([A-Z][a-z]*)([0-9./]*)([+-][0-9]*)?'
    def __init__( self, chemical_formula ):
        if type(chemical_formula) is str:
            self.atoms = self.parse_chemical_formula( chemical_formula )
        else:
            self.atoms = chemical_formula
    def parse_chemical_formula( self, chemical_formula ):
        if type(chemical_formula) is ChemicalFormula:
            chemical_formula = str(chemical_formula)
        pattern = re.compile(self.chemical_formula_RE)
        atoms = []
        charge = 0
        for atomic_symbol, stoichiometry, charge in pattern.findall( chemical_formula ):
            if stoichiometry == '':
                stoichiometry = Fraction(1)
            else:
                stoichiometry = Fraction(stoichiometry)
            if charge == '':
                charge = 0
            elif len(charge) == 1:
                if charge == '+':
                    charge = 1
                else:
                    charge = -1
            elif charge[0] == '+':
                charge = int(charge[1:])
            elif charge[0] == '-':
                charge = int(charge)
            if stoichiometry != 0:
                atoms.append( dict(symbol=atomic_symbol, atom=periodic.element(atomic_symbol), number=stoichiometry ))
        return dict(molecular_formula=atoms, charge=charge)
    
    def get_charge( self ):
        return self.atoms['charge']
    
    def get_chemical_composition( self ):
        composition = {}
        for atom in self.atoms['molecular_formula']:
            if atom['symbol'] in composition:
                composition[atom['symbol']] += atom['number']
            else:
                composition[atom['symbol']] = atom['number']
        composition['Charge'] = self.atoms['charge']
        return pd.Series(composition)
    
    def __contains__( self, item ):
        for atom in self.atoms['molecular_formula']:
            if item == atom:
                return True
        return False
    def get_molecular_mass( self ):
        count = 0
        for atom in self.atoms['molecular_formula']:
            count += atom['atom'].mass*atom['number']
        return count
    
    def get_number_of_molecule_protons( self ):
        count = 0
        for atom in self.atoms['molecular_formula']:
            count += atom['atom'].atomic*atom['number']
        return count
    def get_number_of_atoms( self, element ):
        count = 0
        for atom in self.atoms['molecular_formula']:
            if element == atom['symbol']:
                count += atom['number']
        return count
    def get_number_of_atom_protons( self, element ):
        count = 0
        for atom in self.atoms['molecular_formula']:
            if element == atom['symbol']:
                count += atom['number']*atom['atom'].atomic
        return count
    def get_atom_mass( self, element ):
        count = 0
        for atom in self.atoms['molecular_formula']:
            if element == atom['symbol']:
                count += atom['number']*atom['atom'].mass
        return count
    
    def __add__( self, other):
        return Mixture( str(self )) + Mixture( str( other ))
    
    def __sub__( self, other ):
        return Mixture( str( self )) - Mixture( str( other ))

    def to_latex( self ):
        out = ''
        for atom in self.atoms['molecular_formula']:
            if atom['number'] == Fraction(1):
                out +=  r'\text{%s}' % atom['symbol'] 
            elif atom['number'].denominator == 1:
                out +=  r'\text{%s}_{%d}' % (atom['symbol'], atom['number'].numerator)
            else:
                out +=  r'\text{%s}_\frac{%d}{%d}' % (atom['symbol'], atom['number'].numeratoratom['number'].denominator)
        if self.atoms['charge'] == 0:
            return out
        elif self.atoms['charge'] == -1:
            return out + '^-'
        elif self.atoms['charge'] == 1:
            return out + '^+'
        elif self.atoms['charge'] > 0:
            return out + '^{' + str( self.atoms['charge']) + '+}'
        else:
            return out + '^{' + str(-self.atoms['charge']) + '-}'

    def __str__(self,latex=False):
        out = ''
        for atom in self.atoms['molecular_formula']:
            if atom['number'] == Fraction(1):
                out += atom['symbol']
            elif atom['number'].denominator == 1:
                out += atom['symbol'] + str(atom['number'].numerator)
            elif atom['number'] != 0:
                out += atom['symbol'] + str(atom['number'])
            else:
                out += ''
        if self.atoms['charge'] == 0:
            return out
        elif self.atoms['charge'] == -1:
            return out + '-'
        elif self.atoms['charge'] == 1:
            return out + '+'
        elif self.atoms['charge'] > 1:
            return out + '+' + str( self.atoms['charge'])
        else:
            return out + str(self.atoms['charge'])

    def __repr__( self ):
        return self.__str__()
