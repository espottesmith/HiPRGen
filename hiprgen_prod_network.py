import os
import sys
import subprocess
import sqlite3
import pickle

from monty.serialization import loadfn, dumpfn

from HiPRGen.network_loader import NetworkLoader
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge
from HiPRGen.species_filter import species_filter
from HiPRGen.bucketing import bucket
from HiPRGen.report_generator import ReportGenerator
from HiPRGen.initial_state import insert_initial_state
from HiPRGen.constants import ROOM_TEMP, Terminal
from HiPRGen.reaction_filter_payloads import (
    DispatcherPayload,
    WorkerPayload
)
from HiPRGen.species_questions import *
from HiPRGen.reaction_questions import (
    li_ec_reaction_decision_tree,

)

if len(sys.argv) != 2:
    print("usage: python test.py number_of_threads")
    quit()

number_of_threads = sys.argv[1]

infsep_mols = loadfn("/Users/ewcss/software/lbnl-scripts/analysis/kmc_paper/final/infsep_mols.json")

libe_format = list()
for ii, ifmol in enumerate(infsep_mols):
    this = dict()
    this["molecule_id"] = str(ii)
    this["molecule"] = ifmol["molecule"]
    this["spin_multiplicity"] = this["molecule"].spin_multiplicity
    this["molecule_graph"] = ifmol["molecule_graph"]
    this["thermo"] = {"raw": {"electronic_energy_Ha": ifmol["energy"],
                              "total_enthalpy_kcal/mol": ifmol["enthalpy"],
                              "total_entropy_cal/molK": ifmol["entropy"]}}

    this["partial_charges"] = {"resp": ifmol["resp"]}
    if this["spin_multiplicity"] == 1:
        this["partial_charges"]["mulliken"] = ifmol["mulliken"]
    else:
        this["partial_charges"]["mulliken"] = [m[0] for m in ifmol["mulliken"]]

    libe_format.append(this)



species_decision_tree =  [
    (li_fix_hydrogen_bonding(), Terminal.KEEP),
    (li_set_solvation_free_energy(li_ec), Terminal.KEEP),
    (charge_too_big(), Terminal.DISCARD),
    # (li0_filter(), Terminal.DISCARD),
    (compute_graph_hashes, Terminal.KEEP),
    (metal_ion_filter(), Terminal.DISCARD),
    (bad_metal_coordination(), Terminal.DISCARD),
    (mol_not_connected(), Terminal.DISCARD),
    (metal_complex(), Terminal.DISCARD),
    (add_star_hashes(), Terminal.KEEP),
    (add_unbroken_fragment(), Terminal.KEEP),
    (add_single_bond_fragments(), Terminal.KEEP),
    (species_default_true(), Terminal.KEEP)
    ]

folder = "/Users/ewcss/data/kmc/natchem/hiprgen"

mol_entries = species_filter(
    libe_format,
    mol_entries_pickle_location=folder + '/mol_entries.pickle',
    species_report=folder + '/unfiltered_species_report.tex',
    species_decision_tree=species_decision_tree,
    coordimer_weight=lambda mol: mol.solvation_free_energy
)

bucket(mol_entries, folder + '/buckets.sqlite')

params = {
    'temperature' : ROOM_TEMP,
    'electron_free_energy' : -1.4
}

dispatcher_payload = DispatcherPayload(
    folder + '/buckets.sqlite',
    folder + '/rn.sqlite',
    folder + '/reaction_report.tex'
)

worker_payload = WorkerPayload(
    folder + '/buckets.sqlite',
    li_ec_reaction_decision_tree,
    params,
    Terminal.DISCARD
)

dumpfn(dispatcher_payload, folder + '/dispatcher_payload.json')
dumpfn(worker_payload, folder + '/worker_payload.json')

subprocess.run(
    [
        'mpiexec',
        '--use-hwthread-cpus',
        '-n',
        number_of_threads,
        'python',
        '/Users/ewcss/software/HiPRGen/run_network_generation.py',
        folder + '/mol_entries.pickle',
        folder + '/dispatcher_payload.json',
        folder + '/worker_payload.json'
    ]
)
