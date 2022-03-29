from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def score_similarity(query_smi: str, ref_smi: str, metric=None) -> float:
    query_mol = Chem.MolFromSmiles(query_smi)
    ref_mol = Chem.MolFromSmiles(ref_smi)

    query_fp = AllChem.GetMorganFingerprint(query_mol, 2, useFeatures=True)
    ref_mol = AllChem.GetMorganFingerprint(ref_mol, 2, useFeatures=True)

    score = DataStructs.DiceSimilarity(query_fp, ref_mol)

    return score


def sdf_to_smiles(sdf_data):
    suppl = Chem.SDMolSupplier(sdf_data)
    smiles_list = []
    for mol in suppl:
        smiles = Chem.MolToSmiles(mol)
        smiles_list.append(smiles)

    return smiles_list
