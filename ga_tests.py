import cPickle
import unittest 
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


smiles_string = 'c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1'
correct_coords = [('C', ['-6.0770', '0.7025', '-0.3583']), ('C', ['-5.9729', '-0.3201', '-1.3039']), ('C', ['-4.7242', '-0.8653', '-1.6146']), ('C', ['-3.5664', '-0.3892', '-0.9788']), ('C', ['-2.3086', '-0.9292', '-1.2832']), ('C', ['-1.1540', '-0.4506', '-0.6455']), ('C', ['0.1046', '-0.9904', '-0.9496']), ('C', ['1.2596', '-0.5126', '-0.3126']), ('C', ['2.5180', '-1.0531', '-0.6174']), ('C', ['3.6720', '-0.5751', '0.0197']), ('C', ['4.9332', '-1.1132', '-0.2824']), ('C', ['6.0770', '-0.6296', '0.3583']), ('C', ['5.9729', '0.3930', '1.3039']), ('C', ['4.7242', '0.9381', '1.6146']), ('C', ['3.5664', '0.4621', '0.9788']), ('C', ['2.3086', '1.0021', '1.2832']), ('C', ['1.1540', '0.5235', '0.6455']), ('C', ['-0.1046', '1.0634', '0.9496']), ('C', ['-1.2595', '0.5855', '0.3126']), ('C', ['-2.5180', '1.1261', '0.6174']), ('C', ['-3.6720', '0.6480', '-0.0197']), ('C', ['-4.9332', '1.1861', '0.2824']), ('H', ['-7.0495', '1.1223', '-0.1208']), ('H', ['-6.8647', '-0.6921', '-1.7986']), ('H', ['-4.6706', '-1.6599', '-2.3538']), ('H', ['-2.2287', '-1.7257', '-2.0199']), ('H', ['0.1858', '-1.7870', '-1.6862']), ('H', ['2.6004', '-1.8498', '-1.3538']), ('H', ['5.0417', '-1.9093', '-1.0141']), ('H', ['7.0495', '-1.0494', '0.1208']), ('H', ['6.8647', '0.7650', '1.7986']), ('H', ['4.6706', '1.7328', '2.3538']), ('H', ['2.2287', '1.7986', '2.0199']), ('H', ['-0.1858', '1.8600', '1.6862']), ('H', ['-2.6004', '1.9227', '1.3538']), ('H', ['-5.0417', '1.9822', '1.0141'])]
mol_result_list = cPickle.load( open("100_3_gen_results.p", "rb" ) )



class ga_setup_quick_test(unittest.TestCase):
    

   def test_mutator(self):
       from ga_setup_selection import mutator
       print "testing mutator"
       for i in range(0,10):
           smiles_list = list(smiles_string) #converting to list for mutator
           new_smiles_list = mutator(smiles_list, mutation_factor = 10, N_number = 5)
           new_smiles = ''.join(new_smiles_list)
           new_mol = Chem.MolFromSmiles(new_smiles)
           self.assertIsInstance(new_mol,rdkit.Chem.rdchem.Mol)

   
   def test_coords_gen(self):
       from ga_setup_selection import coords_generator
       print "testing coords gen"
       mol_list = []
       mol = Chem.MolFromSmiles(smiles_string)
       mol_list.append(mol)
       coords = coords_generator(mol_list)[0]
       len1 = len(coords)
       actual_len = len(correct_coords)
       self.assertEqual(len1,actual_len)
        
    def test_tournament(self):
        from ga_setup_selection import tournament_selection
        print "testing tournament selection"         



if __name__ == '__main__':
    unittest.main()
