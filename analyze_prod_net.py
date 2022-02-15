from HiPRGen.network_loader import NetworkLoader
from HiPRGen.report_generator import ReportGenerator

from HiPRGen.mc_analysis import (
    reaction_tally_report,
    species_report,
    Pathfinding,
    SimulationReplayer,
    generate_pathway_report,
    sink_report,
    consumption_report,
    redox_report
)

folder = "/Users/ewcss/data/kmc/natchem/hiprgen"

network_loader = NetworkLoader(
    folder + '/rn.sqlite',
    folder + '/mol_entries.pickle',
    folder + '/li_ec.sqlite'
    )

network_loader.load_trajectories()
network_loader.load_initial_state()

report_generator = ReportGenerator(
    network_loader.mol_entries,
    folder + '/dummy.tex',
    rebuild_mol_pictures=True)

reaction_tally_report(
    network_loader,
    folder + '/li_ec.tex'
)

network_loader = NetworkLoader(
    folder + '/rn.sqlite',
    folder + '/mol_entries.pickle',
    folder + '/li_ec_co2.sqlite'
    )

network_loader.load_trajectories()
network_loader.load_initial_state()

reaction_tally_report(
    network_loader,
    folder + '/li_ec_co2.tex'
)

network_loader = NetworkLoader(
    folder + '/rn.sqlite',
    folder + '/mol_entries.pickle',
    folder + '/li_ec_h2o.sqlite'
    )

network_loader.load_trajectories()
network_loader.load_initial_state()

reaction_tally_report(
    network_loader,
    folder + '/li_ec_h2o.tex'
)

network_loader = NetworkLoader(
    folder + '/rn.sqlite',
    folder + '/mol_entries.pickle',
    folder + '/li_ec_co2_h2o.sqlite'
    )

network_loader.load_trajectories()
network_loader.load_initial_state()

reaction_tally_report(
    network_loader,
    folder + '/li_ec_co2_h2o.tex'
)