"""
reactions.py
 Automate generation of synthetic workflows for the OpenTrons robotics platform

Handles getting reactants and/or products from an input 
"""

import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw

class Reactant(object):
    """
    A chemical reactant object eg. Caffeine

    Attributes:
        name        The chemical name eg. 1,3,7-trimethylpurine-2,6-dione
        SMILES      The SMILES notation eg. CN1C=NC2=C1C(=O)N(C(=O)N2C)C
        location    The physical location of the reactant from the ChemInventory DB eg. RCaH XChem > RCaH XChem > RCaH_1.27_XChem > Freezer 1 > XCHEM-FR1
        comment     The comment added to the ChemInventory DB to assist with finding the reactant eg. C7 (Matrix position of enamine compound in XCHEM-FR1) 
        mol         The rdkit mol object of the reactant
        MW          The exact molecular weight of the reactant
        solubility  The molar solubility of the reactant  

    """ 

    def __init__(self, name, SMILES, location, solubility=None):
        self.name = name
        self.SMILES = SMILES
        self.location = location
        self.comment = comment
        self.mol = Chem.MolFromSmiles(SMILES)
        self.MW = Descriptors.ExactMolWt(self.mol)
        self.solubility = solubility
        
class Reaction(object):
    """
    A chemical reaction

    Attributes:
        reactants       The list of chemical reactant objects - this is used if the chemist would like to predict the product
        reactionsmiles  The reactant part of the reaction represented as SMILES
        reactantmols    The tuple of reactant rdkit mol objects 
        reactantMWs     The list of exact molecular weights of the reactants
        product         The SMILES notation of the product - this is given if the chemist wants to predict the reactants
    """ 
    def __init__(self, reactants=None, product=None):    
        # Reactants tuples of Compound objects
        if not reactants:
            self.reactants = None
        if reactants:
            self.reactants = [reactant for reactant in reactants]
            self.reactionsmiles = str(reactants[0].SMILES + '.' + reactants[1].SMILES)
            self.reactantmols = tuple([reactant.mol for reactant in reactants])
            self.reactantMWs = [reactant.MW for reactant in reactants]       
            
        self.product = product 
        
    
    def getProduct(self):
        """
        This function uses the IBM api to predict the product formed from reacting
        two reactant
        """
        while self.product is None:
            try:
                # IBM API allows a call to be made every 2s with a mximum of 5 per minute
                time.sleep(30)
                response = rxn4chemistry_wrapper.predict_reaction(self.reactionsmiles)
                time.sleep(30)
                results = rxn4chemistry_wrapper.get_predict_reaction_results(response['prediction_id'])
                rxn_smiles = results['response']['payload']['attempts'][0]['smiles']
                self.product = rxn_smiles.split('>>')[-1]          
            except Exception as e:
                print(e)
  
                
    def getReactions(self):
        """
        Use the IBM API to get some possible retrosynthesis routes
        """
        # Create dummy dictionary to create while loop to catch when status is a SUCCESS
        results = {}
        results['status'] = None
        self.reactions = None
        
        while self.reactions is None:  
            try:
                time.sleep(30)
                response = rxn4chemistry_wrapper.predict_automatic_retrosynthesis(product=self.product)
                while results['status'] != 'SUCCESS': 
                    time.sleep(30)
                    results = rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(response['prediction_id'])
                    self.reactions = results
            except Exception as e:
                print(e)
        return self.reactions


    def drawReaction(self):
        try:
            reactionsmiles = "{}>>{}".format(self.reactionsmiles,self.product)
            print(reactionsmiles)
            return Draw.ReactionToImage(reactionsmiles,useSmiles=True)
        except Exception as e: 
            print(e)
    
    def getReactantAmounts(self, product_mass, mol_equivalents=1, product_yield=1):
        """
        mol_equivalents is the scale of reactant two mols needed
        reagents are rdkit mol objects of any extra reagents needed 
        
        Returns estimated masses of reactants needed to yield product mass in mg 
        """
        self.productmass = product_mass
        self.productyield = product_yield
        try:
            self.productmols = self.productmass / Descriptors.ExactMolWt(self.product)
            if len(self.reactantMWs) == 1:
                self.react_1_mass = (self.productmols * self.reactantMWs[0]) / self.productyield 
                self.react_2_mass = 0
            if len(self.reactantMWs) == 2:
                self.react_1_mass = (self.productmols * self.reactantMWs[0]) / self.productyield
                self.react_2_mass = (self.productmols * mol_equivalents * self.reactantMWs[1]) / self.productyield
        except Exception as e: 
            print(e)  
            




