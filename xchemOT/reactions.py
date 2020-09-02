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
import os
from rxn4chemistry import RXN4ChemistryWrapper

# Setup IBM RxN API
api_key=os.environ['IBM_API_KEY'] 
rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=api_key)
rxn4chemistry_wrapper.create_project('Moonshot_amide_synthesis')

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

    def __init__(self, name, SMILES, location, comments, solubility):
        self.name = name
        self.SMILES = SMILES
        self.location = location
        self.comments = comments
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
        productmass     The expected mass (g) of product to be made
    """ 
    def __init__(self, productmass, reactants=None, product=None):    
        # Reactants tuples of Compound objects
        if not reactants:
            self.reactants = self.getReactions()
        if reactants:
            self.reactants = [reactant for reactant in reactants]
            self.reactionsmiles = str(reactants[0].SMILES + '.' + reactants[1].SMILES)
            self.reactantmols = tuple([reactant.mol for reactant in reactants])
            self.reactantMWs = [reactant.MW for reactant in reactants]
            self.reactantsolubility = [reactant.solubility for reactant in reactants]        
        
        # NB need to have option to start with reactants or product
        if not product:    
            self.product = self.getProduct()
            self.productmol = Chem.MolFromSmiles(self.product)
            self.productMW = Descriptors.ExactMolWt(self.productmol)
        if product:
            self.product = product
            self.productmol = Chem.MolFromSmiles(product)
            self.productMW = Descriptors.ExactMolWt(self.productmol)
        
        self.productmass = productmass
        self.productmols = self.productmass / self.productMW
         
    
    def getProduct(self):
        """
        This function uses the IBM api to predict the product formed from reacting
        two reactants
        """
        self.product = None

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
        return self.product 
  
                
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
            self.reactionsmarts = AllChem.ReactionFromSmarts("{}>>{}".format(self.reactionsmiles,self.product),useSmiles=True)
            self.reactionimage = Draw.ReactionToImage(self.reactionsmarts)
            return Draw.ReactionToImage(self.reactionsmarts)
        except Exception as e: 
            print(e)
    
    def getReactantAmounts(self, mol_equivalents=1, product_yield=1):
        """
        mol_equivalents is the mol equivalents of reactant two needed
        reagents are rdkit mol objects of any extra reagents needed 
        
        Returns estimated masses (g) and volumes (ml) of reactants needed to yield product mass in mg 
        """
        
        try:
            if len(self.reactantMWs) == 1:
                # Do calcs for reactant 1 -> required: mass, mols and volume
                self.react_1_mass = self.productmols * self.reactantMWs[0] / product_yield
                self.react_1_mols = self.react_1_mass / self.reactantMWs[0]
                self.react_1_volume = (self.react_1_mols / self.reactantsolubility[0]) * 1000  
                # Do calcs for reactant 2 -> required: mass, mols and volume
                self.react_2_mass = None
                self.react_2_mols = None
                self.react_2_volume = None

            if len(self.reactantMWs) == 2:
                # Do calcs for reactant 1 -> required: mass, mols and volume
                self.react_1_mass = self.productmols * self.reactantMWs[0] / product_yield
                self.react_1_mols = self.react_1_mass / self.reactantMWs[0] 
                self.react_1_volume = (self.react_1_mols / self.reactantsolubility[0]) * 1000  
                # Do calcs for reactant 2 -> required:  mass, mols and volume
                self.react_2_mass = self.react_1_mols * mol_equivalents * self.reactantMWs[1] / product_yield                
                self.react_2_mols = self.react_2_mass / self.reactantMWs[1] 
                self.react_2_volume = (self.react_2_mols / self.reactantsolubility[1]) * 1000 

        except Exception as e: 
            print(e) 

    # Need to write all of this info to Pandas df as some point
    def getDictionary(self):
        return {
            'Product_SMILES': self.product,
            'Product_mass': self.productmass,
            'React_1_name': self.reactants[0].name,
            'React_1_location': self.reactants[0].location + '__' + self.reactants[0].comments,
            'React_1_mass': self.react_1_mass,
            'React_1_vol': self.react_1_volume,
            'React_2_name': self.reactants[1].name,
            'React_2_location': self.reactants[1].location + '__' + self.reactants[1].comments,
            'React_2_mass': self.react_2_mass,
            'React_2_vol': self.react_2_volume          
        } 
            




