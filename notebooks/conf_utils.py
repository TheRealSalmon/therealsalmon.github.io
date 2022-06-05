from typing import Union
from rdkit import Chem
from rdkit.Chem import (
    rdMolDescriptors,
    rdDistGeom,
    AllChem,
)
import py3Dmol


def get_low_energy_conformer(
    input_mol: Union[str, Chem.rdchem.Mol],
    max_iters: int = 200,
    allow_more_rot_bonds: bool = False,
) -> Chem.rdchem.Mol:
    """
    Obtain the lowest energy conformer.

    Finds the MMFF low-energy conformer. It generates n conformers, where n
    depends on the number of rotatable bonds. Then the conformers are optimized
    with the MMFF forcefield. Finally, the lowest energy conformer is returned.
    Will raise error if the number of rotatable bonds is greater than 10.

    Examples
    --------
    mol = Chem.MolFromSmiles('OCCCO')
    low_e_mol = get_low_energy_conformer(mol)

    Parameters
    ----------
    input_mol: `rdkit.Chem.rdchem.Mol`
        The input RDKit mol object.
    max_iters: `int`, default = 200
        The number of iterations allowed for the MMFF optimization.
    allow_more_rot_bonds: bool, default=False
        Flag to allow more than 10 rotatable bonds. Compounds with many
        rotatable bonds have a high chance of failure.

    Returns
    -------
    `rdkit.Chem.rdchem.Mol`
        An RDKit Mol object embedded with the (hopefully) lowest energy
        conformer
    """
    if type(input_mol) == str:
        try:
            mol = Chem.MolFromSmiles(input_mol)
        except ValueError as err:
            print(err)
            raise
    elif type(input_mol) == Chem.rdchem.Mol:
        # create a quickCopy of a Mol which removes all conformers/properties
        mol = Chem.rdchem.Mol(input_mol, True)
    else:
        raise TypeError('expected SMILES string or RDKit Mol')
    mol = Chem.AddHs(mol)
    # we will later return low_e_mol
    low_e_mol = Chem.rdchem.Mol(mol)

    # this count of rotatable bonds INCLUDES amides
    n_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    # this is the latest set of parameters for the ETKDG algorithm
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = 42
    ps.pruneRmsThresh = 0.5
    # we create 100 conformers for 0-4 rotatable bonds, 200 conformers for 5-7
    #  rotatable bonds, 300 conformers for 8-10 rotatablebonds, and if allowed,
    #  500 conformers for >10 rotatable bonds
    if n_rot_bonds <= 4:
        rdDistGeom.EmbedMultipleConfs(mol, 100, ps)
    elif n_rot_bonds > 4 and n_rot_bonds <= 7:
        rdDistGeom.EmbedMultipleConfs(mol, 200, ps)
    elif n_rot_bonds > 7 and n_rot_bonds <= 10:
        rdDistGeom.EmbedMultipleConfs(mol, 300, ps)
    elif n_rot_bonds > 10 and allow_more_rot_bonds:
        rdDistGeom.EmbedMultipleConfs(mol, 500, ps)
    elif n_rot_bonds > 10 and not allow_more_rot_bonds:
        raise ValueError('compounds contains greater than 10 rotatable bonds, '
                         'set allow_more_rot_bonds = True')
    
    tuple_success_energy = AllChem.MMFFOptimizeMoleculeConfs(mol)
    # here we unpack our results
    conformers = []
    energies = []
    for i, tup in enumerate(tuple_success_energy):
        if tup[0] == 0:
            conformers.append(mol.GetConformer(i))
            energies.append(tup[1])
    if len(conformers) == 0:
        raise ValueError('MMFF94 optimization failed, raise max_iters')
    # sorts conformers and energies together on the energies
    energies, conformers = zip(*sorted(zip(energies, conformers)))

    low_e_mol.AddConformer(conformers[0])
    return low_e_mol

def remove_nonpolar_hs(input_mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    """Remove nonpolar hydrogen atoms.

    Finds all hydrogens bonded to carbon atoms and returns an RDKit Mol object
    with these hydrogens removed.

    Examples
    --------
    mol = Chem.MolFromSmiles('OCCCO')
    mol_polar_h = remove_nonpolar_hs(mol)

    Parameters
    ----------
    input_mol: `rdkit.Chem.rdchem.Mol`
        The input RDKit mol object.

    Returns
    -------
    `rdkit.Chem.rdchem.Mol`
        An RDKit Mol object with all nonpolar hydrogens removed."""
    # Make a copy of input mol.
    mol = Chem.rdchem.Mol(input_mol)

    # Find indices of all hydrogens bonded to carbons.
    nonpolar_hs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            for n in atom.GetNeighbors():
                if n.GetAtomicNum() == 1:
                    nonpolar_hs.append(n.GetIdx())
    # The list needs to be ordered from high-to-low to avoid indexing issues.
    nonpolar_hs = sorted(nonpolar_hs, reverse=True)

    # We create a Read/Write Mol and remove the nonpolar hydrogens.
    rwmol = Chem.rdchem.RWMol(mol)
    for h in nonpolar_hs:
        rwmol.RemoveAtom(h)

    return rwmol.GetMol()

def display_3d_mol(
    mol: Chem.rdchem.Mol,
    nonpolar_h: bool = False,
    carbon_color: str = 'random',
) -> None:
    """Use py3Dmol to visualize mol in 3D.

    Examples
    --------
    mol = Chem.MolFromSmiles('OCCCO')
    mol = get_low_energy_conformer(mol)
    display_3d_mol(mol)

    Parameters
    ----------
    mol: `rdkit.Chem.rdchem.Mol`
        The input RDKit mol object with an embedded 3D conformer.
    nonpolar_h: `bool`, default = False
        Whether or not to show nonpolar (C-H) hydrogens"""
    from random import choice
    
    mol_block = ''
    if nonpolar_h:
        mol_block = Chem.rdmolfiles.MolToMolBlock(mol, includeStereo=True)
    else:
        mol_block = Chem.rdmolfiles.MolToMolBlock(remove_nonpolar_hs(mol),
                                                  includeStereo=True)
    if carbon_color == 'random':
        carbon_color = choice([
            'white', 'pink', 'cyan', 'green', 'magenta', 'yellow', 'orange'
        ])
    elif carbon_color in [
        'white', 'pink', 'cyan', 'green', 'magenta', 'yellow', 'orange'
    ]:
        pass
    else:
        raise ValueError(f'{carbon_color} is not a valid color')
    view = py3Dmol.view(data=mol_block,
                        width=400,
                        height=300,
                        style={
                            'stick': {
                                'colorscheme': f'{carbon_color}Carbon',
                                'radius':0.25
                            }
                        })
    view.setViewStyle({'style':'outline','color':carbon_color,'width':0.04})
    view.setBackgroundColor('#111111')
    view.show()
