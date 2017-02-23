# RxnEqn
This is code to support the computation of stoichiometric equations for the Thermodynamic Electon Equivalents Model (TEEM) developed by Rittman and McCarty.

## Installation
To install:

```
git clone https://github.com/djinnome/rxneqn.git
cd rxneqn
pip install .
```

## Usage


### ChemicalFormula 
```python
chemical_formula = 'C6H12O5.5N'
cf = ChemicalFormula(chemical_formula)
cf.get_chemical_composition()
```

    C            6
    Charge       0
    H           12
    N            1
    O         11/2
    dtype: object

### Mixture

```python
mix = Mixture('3 C6H12O6 + 4 NO2-')
mix.subtract_from_mixture(Mixture('5 NO2- + H+'))
```


    3 C6H12O6

### Reaction

```python
rxn = Reaction('1/4 CO2 + 1/12NH3 + H+  ==> 1/12 CH3CHNH2COOH + 1/3 H2O + OH-')
rxn -  Reaction('CO2 + H2O ==> 4 HCO3')
```
    H+ + 4 HCO3 + 1/12 NH3 ==> 1/12 CH3CHNH2COOH + 3/4 CO2 + 4/3 H2O + OH-


```python
rxn.normalize()
```

    3 CO2 + 12 H+ + NH3 ==> CH3CHNH2COOH + 4 H2O + 12 OH-

### HalfReactionBalancer

```python
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
1/12 CO2 + E- + H+ + 1/6 HCO3- + 1/12 NH4+ ==> 1/12 CH3CHNH2COO- + 1/2 H2O""".split('\n')
eqn = {}
eqn[13] = half_rxn.balance_half_reaction('C6H12O6', 'CO2')
eqn[14] = half_rxn.balance_half_reaction('NO3-', 'N2',nitrogen_source='NO3-')
eqn[15] = eqn[13] +  eqn[14] 
eqn[16] = half_rxn.balance_half_reaction(carbon_dioxide, alanine, nitrogen_source)
eqn[17] = half_rxn.balance_half_reaction('CO2', 'C5H7O2N','NH4+')
table3 = {'O-12': '1/8 CO2 + H+ + E- ==> 1/8 CH4 + 1/4 H2O',
         'O-2': '1/6 CO2 + 1/12 NH4+ + 1/12 HCO3- + H+ + E- ==> 1/12 CH3CHNH2COOH + 5/12 H2O',
         'O-3': str(half_rxn.balance_half_reaction('CO2','C6H5COO-'))}

eqn[18] = Reaction(table3['O-12'])
eqn[19] = Reaction(table3['O-2'])
eqn[20] = eqn[18] - eqn[19]
eqn[21] = half_rxn.balance_half_reaction('CO2','CH3COCOO-')
eqn[22] = half_rxn.balance_half_reaction('CO2','CH3CHNH2COOH','NH3')
eqn[23] = half_rxn.balance_half_reaction('CO2','CH3CHNH2COOH','NH4+')
eqn[24] = half_rxn.balance_half_reaction('CO2','CH3CHNH2COO-','NH4+')

k = 0

for i in range(13,25):
    assert str(eqn[i]) == truth_eqn[k]
    k += 1
print("{} tests passed".format(k))
```

    12 tests passed


