"""Regenerate the figures of docs/source/tutorials/visualization.rst (run from the repo root)."""

import os
import sys
import warnings

sys.path.insert(0, "src")
warnings.simplefilter("ignore")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from phylozoo import DirectedPhyNetwork  # noqa: E402
from phylozoo.core.network.dnetwork.derivations import to_sd_network  # noqa: E402
from phylozoo.viz import plot  # noqa: E402
from phylozoo.viz.sdnetwork.style import SDNetStyle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
NEWICK = (
    "((Ps.jonesii,((((X.evelynae,(X.xiphidium,(X.maculatus)#H1)),"
    "(((X.meyeri,X.couchianus),X.gordoni),X.milleri)),X.variatus),"
    "((#H1,((((((X.multilineatus,X.nigrensis),(X.pygmaeus,X.continens)),"
    "((X.malinche,X.birchmanni),(X.cortezi)#H2)),(X.montezumae,"
    "(#H2,X.nezahualcoyotl))),((X.signum,(X.hellerii,(X.mayae,"
    "X.alvarezi))))#H3),(#H3,(X.clemenciae,X.monticolus)))),X.andersi))));"
)
dnet = DirectedPhyNetwork.from_string(NEWICK, format="enewick")
sd = to_sd_network(dnet)


def save(net, name, size=(9, 7), **kwargs):
    fig, ax = plt.subplots(figsize=size)
    plot(net, ax=ax, **kwargs)
    fig.savefig(os.path.join(HERE, name), dpi=130, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


save(dnet, "d_cladogram.png", size=(11, 6))
save(dnet, "d_cladogram_lr.png", size=(8, 9), direction="LR")
save(dnet, "d_layered_aligned.png", size=(11, 6), layout="pz-layered", align_leaves=True)
save(dnet, "d_radial.png", layout="pz-radial")
save(sd, "sd_unrooted.png")
style = SDNetStyle(
    node_size=60, leaf_size=140, node_color="#dddddd", leaf_color="#1f4e79",
    hybrid_color="#f5b041", edge_color="#555555", hybrid_edge_color="#d35400",
    edge_width=1.6, label_offset=0.012, label_font_size=11,
)
fig, ax = plt.subplots(figsize=(10, 7))
plot(sd, ax=ax, style=style)
fig.savefig(os.path.join(HERE, "sd_styled.png"), dpi=150, bbox_inches="tight", pad_inches=0.05)
print("figures written to", HERE)
print(dnet.number_of_nodes(), "nodes,", len(dnet.taxa), "taxa,", len(dnet.hybrid_nodes), "hybrid nodes")
print(sd.number_of_nodes(), "nodes,", len(sd.taxa), "taxa,", len(sd.hybrid_nodes), "hybrid nodes (sd)")
print(dnet.to_preview_string(max_width=100))
