# Copyright 2026 EMBL - European Bioinformatics Institute
# Author: Sri Devan Appasamy
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Interactive residue-pair explorer (ipywidgets).

One entry point, `residue_pair_explorer`, which displays a filter control and
renders the residue-pair table and frequency heatmap for the selected scope.
Requires a live kernel; in a non-interactive execution the widgets render as
static output.
"""

from __future__ import annotations

from collections import Counter

import ipywidgets as widgets
from IPython.display import clear_output, display

from pdbe_interfaces import outputs, plots
from pdbe_interfaces.representation import InterfaceRecord
from pdbe_interfaces.similarity import ClusterResult


def residue_pair_explorer(
    records: list[InterfaceRecord],
    cluster_result: ClusterResult,
    partner_map: dict[tuple[str, int], str] | None = None,
    top_n: int = 15,
) -> None:
    """Display the residue-pair table and heatmap with a cluster filter.

    The dropdown selects the whole dataset or a single interface interaction
    state; the slider caps the table rows and the heatmap dimension on each
    side. Selecting one state shows the residue-pair composition of that state
    alone.
    """
    assignment = cluster_result.flat_assignment
    sizes = Counter(assignment.tolist())

    options = [("All clusters", "all")]
    options += [
        (f"Cluster {cid} (n={size})", str(cid))
        for cid, size in sorted(sizes.items(), key=lambda kv: (-kv[1], kv[0]))
    ]

    cluster_dropdown = widgets.Dropdown(
        options=options, value="all", description="Filter:",
        style={"description_width": "initial"},
    )
    top_n_slider = widgets.IntSlider(
        value=top_n, min=5, max=40, step=5,
        description="Top N rows / heatmap size",
        style={"description_width": "initial"},
        layout=widgets.Layout(width="500px"),
    )
    table_out = widgets.Output(layout=widgets.Layout(width="42%", overflow="auto"))
    plot_out = widgets.Output(layout=widgets.Layout(width="58%", overflow="auto"))

    def _subset() -> tuple[list[InterfaceRecord], str]:
        selection = cluster_dropdown.value
        if selection == "all":
            return records, f"all {len(records)} interfaces"
        cid = int(selection)
        members = [r for i, r in enumerate(records) if assignment[i] == cid]
        return members, f"cluster {cid} only ({len(members)} interfaces)"

    def _render(_change=None) -> None:
        subset, scope = _subset()
        n = top_n_slider.value
        freq = (
            outputs.interface_frequency_summary(subset, partner_map=partner_map)
            if subset else None
        )
        with table_out:
            clear_output(wait=True)
            if freq is None:
                print("No interfaces in this selection.")
            else:
                print(f"Top {n} residue pairs ({scope})")
                display(freq["pairs"][["contact", "n_interfaces", "fraction"]].head(n))
        with plot_out:
            clear_output(wait=True)
            if freq is not None:
                plots.pair_frequency_heatmap(freq, n, scope)

    cluster_dropdown.observe(_render, names="value")
    top_n_slider.observe(_render, names="value")
    display(widgets.VBox([
        widgets.HBox([cluster_dropdown, top_n_slider]),
        widgets.HBox([table_out, plot_out]),
    ]))
    _render()
