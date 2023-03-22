from rdkit import Chem


def convert_to_mol(smiles):
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol:
        molblock = Chem.MolToMolBlock(rdmol)
        return molblock

    return None


def get_inchi_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        inchi = Chem.MolToInchi(mol)
        inchi_key = Chem.InchiToInchiKey(inchi)
        return inchi_key

    return None


def get_truncated_inchi_key(inchi_key):
    if inchi_key:
        return inchi_key[0:14]
    
    return None


def get_standard_smiles(standardizer_output):
    mol = Chem.MolFromMolBlock(standardizer_output)
    smiles = Chem.MolToSmiles(mol)
    return smiles
