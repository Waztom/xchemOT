"""
Unit and regression test for the xchemOT package.
"""

# Import package, test suite, and other packages as needed
import xchemOT
from xchemOT.reactions import Reactant
import pytest
import sys

def test_xchemOT_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "xchemOT" in sys.modules
    
def test_reactant_class():
    reactant = Reactant(name='1,3,7-trimethylpurine-2,6-dione', 
                        SMILES='CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 
                        location='coffee cup',
                        comments='wake me up',        
                        solubility='100M')

    assert reactant.MW == 194.08037556 