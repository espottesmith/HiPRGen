from HiPRGen.mol_entry import MoleculeEntry, FragmentComplex
import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import copy
from functools import partial
from HiPRGen.constants import li_ec, Terminal, mg_g2, mg_thf, m_formulas, metals
import numpy as np
from monty.json import MSONable

"""
species decision tree:

A question is a function q(mol_entry) -> Bool

Unlike for reaction filtering, these questions should not modify the mol_entry in any way.

A node is either a Terminal or a non empty list [(question, node)]

class Terminal(Enum):
    KEEP = 1
    DISCARD = -1

For the return value of a question, True means travel to this node and False means try next question in the list.

for non terminal nodes, it is an error if every question returns False. i.e getting stuck at a non terminal node is an error.

Once a Terminal node is reached, it tells us whether to keep or discard the species.
"""

def run_decision_tree(mol_entry,
                      decision_tree,
                      decision_pathway=None):

    node = decision_tree

    while type(node) == list:
        next_node = None
        for (question, new_node) in node:
            if question(mol_entry):

                # if decision_pathway is a list,
                # append the question which
                # answered true i.e the edge we follow
                if decision_pathway is not None:
                    decision_pathway.append(question)


                next_node = new_node
                break

        node = next_node


    if type(node) == Terminal:
        if decision_pathway is not None:
            decision_pathway.append(node)


        if node == Terminal.KEEP:
            return True
        else:
            return False
    else:
        print(node)
        raise Exception("unexpected node type reached")


class metal_ion_filter(MSONable):
    "only allow positively charged metal ions"
    def __init__(self):
        pass

    def __call__(self, mol_entry):
        if mol_entry.formula in m_formulas and mol_entry.charge <= 0:
            return True
        else:
            return False

class mol_not_connected(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        return not nx.is_connected(mol.graph)


class add_star_hashes(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        for i in range(mol.num_atoms):
            if i not in mol.m_inds:
                neighborhood = nx.generators.ego.ego_graph(
                    mol.covalent_graph,
                    i,
                    1,
                    undirected=True)

                mol.star_hashes[i] = weisfeiler_lehman_graph_hash(
                    neighborhood,
                    node_attr='specie')

        return False

class add_unbroken_fragment(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        if mol.formula in m_formulas:
            return False

        fragment_complex = FragmentComplex(
             1,
             0,
             [],
             [mol.covalent_hash])

        mol.fragment_data.append(fragment_complex)

        return False

class add_single_bond_fragments(MSONable):

    def __init__(self):
        pass

    def __call__(self, mol):

        if mol.formula in m_formulas:
            return False



        for edge in mol.covalent_graph.edges:
            fragments = []
            h = copy.deepcopy(mol.covalent_graph)
            h.remove_edge(*edge)
            connected_components = nx.algorithms.components.connected_components(h)
            for c in connected_components:

                subgraph = h.subgraph(c)

                fragment_hash = weisfeiler_lehman_graph_hash(
                    subgraph,
                    node_attr='specie')


                fragments.append(fragment_hash)

            fragment_complex = FragmentComplex(
                len(fragments),
                1,
                [edge[0:2]],
                fragments)

            mol.fragment_data.append(fragment_complex)

        return False


class metal_complex(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        # if mol is a metal, it isn't a metal complex
        if mol.formula in m_formulas:
            return False

        return not nx.is_connected(mol.covalent_graph)

class li_fix_hydrogen_bonding(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        if mol.num_atoms > 1:
            for i in range(mol.num_atoms):
                if mol.species[i] == 'H':

                    adjacent_atoms = []

                    for bond in mol.graph.edges:
                        if i in bond[0:2]:

                            if i == bond[0]:
                                adjacent_atom = bond[1]
                            else:
                                adjacent_atom = bond[0]

                            displacement = (mol.atom_locations[adjacent_atom] -
                                            mol.atom_locations[i])

                            dist = np.inner(displacement, displacement)

                            adjacent_atoms.append((adjacent_atom, dist))


                    closest_atom, _ = min(adjacent_atoms, key=lambda pair: pair[1])

                    for adjacent_atom, _ in adjacent_atoms:
                        if adjacent_atom != closest_atom:
                            mol.graph.remove_edge(i, adjacent_atom)
                            mol.covalent_graph.remove_edge(i, adjacent_atom)

        return False

class mg_fix_hydrogen_bonding(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        # TODO: Make this user-defined?
        max_dist = 1.5

        if mol.num_atoms > 1:
            for i in range(mol.num_atoms):
                if mol.species[i] == 'H':

                    adjacent_atoms = []

                    for bond in mol.graph.edges:
                        if i in bond[0:2]:

                            if i == bond[0]:
                                adjacent_atom = bond[1]
                            else:
                                adjacent_atom = bond[0]

                            displacement = (mol.atom_locations[adjacent_atom] -
                                            mol.atom_locations[i])

                            dist = np.inner(displacement, displacement)

                            adjacent_atoms.append((adjacent_atom, dist))

                    for adjacent_atom, dist in adjacent_atoms:
                        if dist > max_dist ** 2:
                            mol.graph.remove_edge(i, adjacent_atom)
                            if adjacent_atom in mol.covalent_graph:
                                mol.covalent_graph.remove_edge(i, adjacent_atom)

        return False

class bad_metal_coordination(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):

        if mol.formula not in m_formulas:

            if (len(metals.intersection(set(mol.species))) > 0 and
                mol.number_of_coordination_bonds == 0):

                return True

        return False

class li_set_solvation_free_energy(MSONable):
    """
    metal atoms coordinate with the surrounding solvent. We need to correct
    free energy to take this into account. The correction is
    solvation_correction * (
           max_coodination_bonds -
           number_of_coordination_bonds_in_mol).
    Since coordination bonding can't reliably be detected from the molecule
    graph, we search for all atoms within a radius of the metal atom and
    discard them if they are positively charged.
    """

    def __init__(self, solvation_env):
        self.solvation_env = solvation_env

    def __call__(self, mol):

        mol.number_of_coordination_bonds = 0

        correction = 0.0

        for i in mol.m_inds:

            species = mol.species[i]
            coordination_partners = []
            radius = self.solvation_env["coordination_radius"][species]

            for j in range(mol.num_atoms):
                if j != i:
                    displacement_vector = (
                        mol.atom_locations[j] -
                        mol.atom_locations[i])

                    if (np.inner(displacement_vector, displacement_vector)
                        < radius ** 2 and (
                            mol.partial_charges_resp[j] < 0 or
                            mol.partial_charges_mulliken[j] < 0)):
                        if not mol.graph.has_edge(i,j):
                            mol.graph.add_edge(i,j)
                        coordination_partners.append(j)


            number_of_coordination_bonds = len(coordination_partners)
            mol.number_of_coordination_bonds += number_of_coordination_bonds
            correction += self.solvation_env["solvation_correction"][species] * (
                self.solvation_env["max_number_of_coordination_bonds"][species] -
                number_of_coordination_bonds)

        mol.solvation_free_energy = correction + mol.free_energy
        return False


class mg_set_solvation_free_energy(MSONable):
    """
    metal atoms coordinate with the surrounding solvent. We need to correct
    free energy to take this into account. The correction is
    solvation_correction * (
           max_coodination_bonds -
           number_of_coordination_bonds_in_mol).
    Since coordination bonding can't reliably be detected from the molecule
    graph, we search for all atoms within a radius of the metal atom and
    discard them if they are positively charged.
    """

    def __init__(self, solvation_env):
        self.solvation_env = solvation_env

    def __call__(self, mol):
        correction = 0.0
        mol.number_of_coordination_bonds = 0

        for i in mol.m_inds:

            species = mol.species[i]
            partial_charge = max( mol.partial_charges_mulliken[i],
                                  mol.partial_charges_resp[i])

            if partial_charge < 1.2:
                effective_charge = "_1"
            elif partial_charge >= 1.2:
                effective_charge = "_2"


            coordination_partners = list()
            species_charge = species + effective_charge
            radius = self.solvation_env["coordination_radius"][species_charge]

            for j in range(mol.num_atoms):
                if j != i:
                    displacement_vector = (
                        mol.atom_locations[j] -
                        mol.atom_locations[i])
                    if (np.inner(displacement_vector, displacement_vector)
                        < radius ** 2 and (
                            mol.partial_charges_resp[j] < 0 or
                            mol.partial_charges_mulliken[j] < 0)):
                        coordination_partners.append(j)

            number_of_coordination_bonds = len(coordination_partners)
            mol.number_of_coordination_bonds += number_of_coordination_bonds
            correction += self.solvation_env[
                "solvation_correction"][species_charge] * (
                self.solvation_env[
                    "max_number_of_coordination_bonds"][species_charge] -
                number_of_coordination_bonds)

        mol.solvation_free_energy =  correction + mol.free_energy
        return False

class no_bare_mg(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        if mol.formula in m_formulas:
            return True
        else:
            return False

class default_true(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        return True


def compute_graph_hashes(mol):
    mol.total_hash = weisfeiler_lehman_graph_hash(
        mol.graph,
        node_attr='specie')

    mol.covalent_hash = weisfeiler_lehman_graph_hash(
        mol.covalent_graph,
        node_attr='specie')

    return False


class li0_filter(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        # some molecules don't have NBO data
        if not mol.partial_charges_nbo:
            return False

        for i in mol.m_inds:
            if (mol.species[i] == 'Li' and
                mol.partial_charges_nbo[i] < 0.1):
                return True

        return False

class charge_too_big(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        if mol.charge > 1 or mol.charge < -1:
            return True

        else:
            return False

# any species filter which modifies bonding has to come before
# any filter checking for connectivity (which includes the metal-centric complex filter)

li_ec_species_decision_tree = [
    (li_fix_hydrogen_bonding(), Terminal.KEEP),
    (li_set_solvation_free_energy(li_ec), Terminal.KEEP),
    (charge_too_big(), Terminal.DISCARD),
    (li0_filter(), Terminal.DISCARD),
    (compute_graph_hashes, Terminal.KEEP),
    (metal_ion_filter(), Terminal.DISCARD),
    (bad_metal_coordination(), Terminal.DISCARD),
    (mol_not_connected(), Terminal.DISCARD),
    (metal_complex(), Terminal.DISCARD),
    (add_star_hashes(), Terminal.KEEP),
    (add_unbroken_fragment(), Terminal.KEEP),
    (add_single_bond_fragments(), Terminal.KEEP),
    (default_true(), Terminal.KEEP)
    ]

mg_g2_species_decision_tree = [
    (mg_set_solvation_free_energy(mg_g2), Terminal.KEEP),
    (no_bare_mg(), Terminal.DISCARD),
    (mg_fix_hydrogen_bonding(), Terminal.KEEP),
    (compute_graph_hashes, Terminal.KEEP),
    (metal_ion_filter(), Terminal.DISCARD),
    (bad_metal_coordination(), Terminal.DISCARD),
    (mol_not_connected(), Terminal.DISCARD),
    (metal_complex(), Terminal.DISCARD),
    (add_star_hashes(), Terminal.KEEP),
    (add_unbroken_fragment(), Terminal.KEEP),
    (add_single_bond_fragments(), Terminal.KEEP),
    (default_true(), Terminal.KEEP)
    ]

mg_thf_species_decision_tree = [
    (mg_set_solvation_free_energy(mg_thf), Terminal.KEEP),
    (no_bare_mg(), Terminal.DISCARD),
    (mg_fix_hydrogen_bonding(), Terminal.KEEP),
    (compute_graph_hashes, Terminal.KEEP),
    (metal_ion_filter(), Terminal.DISCARD),
    (bad_metal_coordination(), Terminal.DISCARD),
    (mol_not_connected(), Terminal.DISCARD),
    (metal_complex(), Terminal.DISCARD),
    (add_star_hashes(), Terminal.KEEP),
    (add_unbroken_fragment(), Terminal.KEEP),
    (add_single_bond_fragments(), Terminal.KEEP),
    (default_true(), Terminal.KEEP)
    ]
