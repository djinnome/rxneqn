from . import Reaction, Mixture
from .utils import LCM,  LCD
import pandas as pd
from fractions import Fraction
class HalfReactionBalancer:
    def __init__( self ):
        pass

    def custom_half_reaction( self, C, H, O, N, charge=0):
        """generate custom half reaction from empirical formula
        (n-c)/d CO2 + c/d NH4+ + c/d HCO3- + (d+f)/d H+ + E- ==> 1/d CnHaObNc  + (2*n -b + c)/d H2O
where d = (4*n + a  - 2*b  - 3*c)
        """
        n,a,b,c, f = Fraction(C),Fraction(H),Fraction(O),Fraction(N), charge
        d = (4*n + a - 2*b -3*c)
        if f == 0:
            biomass_charge = ''
        elif f == 1:
            biomass_charge = '+'
        elif f > 1:
            biomass_charge = '+{}'.format(f)
        else:
            biomass_charge = '{}'.format(f)
        stoichiometry = dict(n=C,a=H, b=O,c=N,
                             CO2=(n-c)/d, NH4 = c/d, HCO3 = c/d,
                             biomass=1/d,H2O = (2*n - b + c)/d, proton=(d+f)/d,
                             charge=biomass_charge)
        eqn = '{CO2} CO2 + {NH4} NH4+ + {HCO3} HCO3- + {proton} H+ + E- ==> {biomass} C{n}H{a}O{b}N{c}{charge} + {H2O} H2O'
        return Reaction(eqn.format(**stoichiometry))
    def balance_half_reaction( self, oxidized_form, reduced_form, nitrogen_source='NH3' ):
        return self.normalize_by_electron(
                    self.balance_charge(
                        self.balance_hydrogen(
                            self.balance_oxygen(
                                self.balance_nonwater_atoms(
                                    self.add_species( 
                                        self.setup_reduction(
                                            oxidized_form, reduced_form),
                                        nitrogen_source))))))
    def setup_reduction( self, oxidized_form, reduced_form ):
        return Reaction(str(oxidized_form) + ' ==> ' + str(reduced_form) )
    
    def balance_element( self, rxn, element_to_be_balanced ):
        rxn = Reaction( str(rxn))
        molecules_of_element = rxn.get_chemical_composition().loc[element_to_be_balanced].fillna(0)
        molecules_of_element = molecules_of_element[molecules_of_element!=0]
        lcm = LCM(molecules_of_element)
        stoich = lcm/molecules_of_element
        return rxn.multiply_factor( stoich )

    def add_species( self, rxn1, nitrogen_source ):
        if 'N' in rxn1.get_chemical_composition().index:
            reactant2 = rxn1.rxn['reactant'].add_to_mixture(Mixture('H2O + ' + str(nitrogen_source)))
        else:
            reactant2 = rxn1.rxn['reactant'].add_to_mixture(Mixture('H2O'))
        product2 = Mixture(str(rxn1.rxn['product']))
        return Reaction( str(reactant2) + ' ==> ' + str(product2))
    
    def balance_nonwater_atoms( self, rxn2 ):
        step3 = Reaction( str( rxn2 ))
        for element in step3.get_chemical_composition().index:
            if element not in ['H', 'O', 'Charge']:
                step3 = self.balance_element( step3, element )
        charge = step3.get_charge()
        if 'C' in step3.get_chemical_composition().index:
            if step3.get_charge() < 0:
                step3.rxn['reactant'] = step3.rxn['reactant'].\
                    add_to_mixture(Mixture('{c:} HCO3-'.format(c=-step3.get_charge()))).\
                    subtract_from_mixture(Mixture('{c:} CO2'.format(c=-step3.get_charge())))
            elif step3.get_charge() > 0:
                step3.rxn['reactant'] = step3.rxn['reactant'].\
                    subtract_from_mixture(Mixture('{c:} HCO3-'.format(c=step3.get_charge()))).\
                    add_to_mixture(Mixture('{c:} CO2'.format(c=step3.get_charge())))

        return step3    
    
    def balance_oxygen( self, rxn3 ):
        step4 = Reaction( str( rxn3 ))
        num_O = step4.rxn['reactant'].get_number_of_atoms('O') - step4.rxn['product'].get_number_of_atoms('O')
        water = '{} H2O ==> {} H2O'
        if num_O > 0:
            return step4.subtract_reaction( Reaction(water.format(num_O, 0)))
        else:
            return step4.subtract_reaction( Reaction(water.format(0, -num_O)))
    
    def balance_hydrogen( self, rxn4 ):
        step5 = Reaction( str( rxn4 ))
        num_H = step5.rxn['reactant'].get_number_of_atoms('H') - step5.rxn['product'].get_number_of_atoms('H')
        protons = '{} H+ ==> {} H+'
        if num_H > 0:
            return step5.subtract_reaction( Reaction( protons.format(num_H, 0)))
        else:
            return step5.subtract_reaction( Reaction( protons.format(0, -num_H)))
    def balance_charge( self, rxn5 ):
        step6 = Reaction( str( rxn5 ))
        
        num_charge = step6.rxn['reactant'].get_charge_of_mixture() - step6.rxn['product'].get_charge_of_mixture()
        charge = '{} E- ==> {} E-'
        if num_charge >0:
            return step6.add_reaction( Reaction( charge.format( num_charge, 0)))
        else:
            return step6.add_reaction( Reaction( charge.format( 0, -num_charge)))
    
    def normalize_by_electron( self, rxn ):
        step7 = Reaction( str( rxn ))
        electrons = rxn.get_stoichiometry_of_species( 'E-')
        factor = pd.Series(Fraction(1,electrons), index=[str(s) for s in step7.get_species()], dtype='object')
        if electrons > 0:
            return step7.multiply_factor( factor )
        elif electrons == 0:
            return step7
        else:
            return step7.multiply_factor( -factor )
