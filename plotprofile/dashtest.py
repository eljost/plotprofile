#!/usr/bin/env python3

from collections import namedtuple
import itertools as it
from pprint import pprint

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import pandas as pd


AU2KJMOL = 2625.499
CAL2J = 4.1868


Rx = namedtuple("Rx", "scf scf_alt dG_solv")


def get_app(energies_df):
    # first_row = energies_df.iloc[0]
    # energies_energies_df.single_point -= first_row.single_point
    # energies_df["single_point"] -= first_row.single_point
    # energies_df["single_point_alt"] -= first_row.single_point_alt
    rx_names = energies_df["name"].unique()
    unique_Ts = energies_df["T"].unique()
    print("unique Ts", unique_Ts)

    rx_radio_options = [{"label": lbl, "value": lbl} for lbl in rx_names]
    temp_slider_kwargs = {
        "min": unique_Ts.min(),
        "max": unique_Ts.max(),
        "step": unique_Ts[1] - unique_Ts[0],
        "value": unique_Ts.min(),
    }

    app = dash.Dash(__name__)
    app.layout = html.Div([
        html.Div([
            html.H2("Energies"),
            dcc.Graph(id="rx_graph"),
            dcc.Graph(id="all_rx_graph"),
            ]),
        html.Div([
            html.H2("Options"),
            html.H3("Sub-reaction"),
            dcc.RadioItems(
                id="rx_radio",
                options=rx_radio_options,
                value=rx_names[0],
            ),
            html.H3("Temperature"),
            dcc.Slider(id="temp_slider", **temp_slider_kwargs),
            html.Div(id="cur_temp"),
            html.H3("Other corrections"),
            dcc.Checklist(
                id="rx_options",
                options=[
                    {"label": "Solvation", "value": "solv"},
                    {"label": "Alt. single point", "value": "alt_sp"},
                ],
                value=["solv", "alt_sp"],
            ),
            html.H3(id="rx_barrier"),
        ], )
    ], style={"columnCount": 2})


    def get_energies(rx_options, rx_name, cur_temp):
        rx_data = energies_df[energies_df.name == rx_name]
        rx_data = rx_data[rx_data["T"] == cur_temp]
        scf = rx_data["single_point"]
        scf_alt = rx_data["single_point_alt"]
        dG_solv = rx_data["dG_solv"]
        dG = rx_data["dG"]

        rx_ens = scf_alt if "alt_sp" in rx_options else scf
        rx_ens += dG
        if "solv" in rx_options:
            rx_ens += dG_solv
        return (rx_ens * AU2KJMOL).values

    def get_figure(rx_options, cur_temp, rx_name=None):
        this_rx_names = [rx_name, ]
        if rx_name is None:
            this_rx_names = rx_names

        # first_row = energies_df.iloc[0]
        # energies_energies_df.single_point -= first_row.single_point
        # energies_df["single_point"] -= first_row.single_point
        # energies_df["single_point_alt"] -= first_row.single_point_alt

        rx_ens = [get_energies(rx_options, name, cur_temp) for name in this_rx_names]
        rx_ens = np.array(rx_ens).flatten()
        rx_ens -= rx_ens[0]

        data = [
            {
                "x": np.arange(rx_ens.size),
                "y": rx_ens,
                "mode": "lines+markers",
                "name": rx_name,
            },
        ]

        layout = {
            # "width": 1100,
            # "height": 450,
            "xaxis": {
                "title": "Rx. coordinate",
            },
            "yaxis": {
                "title": "Energy / kJ mol⁻¹",
                # "range": []
            },
        }

        figure = {
            "data": data,
            "layout": layout,
        }
        return figure

    @app.callback(
        Output("rx_graph", "figure"),
        [Input("rx_options", "value"),
         Input("rx_radio", "value"),
         Input("temp_slider", "value"),]
    )
    def draw_reaction(rx_options, rx_name, cur_temp):
        return get_figure(rx_options, cur_temp, rx_name)

    @app.callback(
        Output("all_rx_graph", "figure"),
        [Input("rx_options", "value"),
         Input("temp_slider", "value"),]
    )
    def draw_all_reactions(rx_options, cur_temp):
        figure = get_figure(rx_options, cur_temp)
        return figure

    @app.callback(
        Output("cur_temp", "children"),
        [Input("temp_slider", "value"),]
    )
    def show_temperature(cur_temp):
        return f"Current temperature: {cur_temp:.2f} K."

    @app.callback(
        Output("rx_barrier", "children"),
        [Input("rx_options", "value"),
         Input("rx_radio", "value"),
         Input("temp_slider", "value"),]
    )
    def set_rx_barrier(rx_options, rx_name, cur_temp):
        rx_ens = get_energies(rx_options, rx_name, cur_temp)
        barrier = rx_ens[1] - rx_ens[0]
        barrier_cal = barrier / CAL2J
        return f"Barrier: {barrier:.1f} kJ mol⁻¹, {barrier_cal:.1f} kcal mol⁻¹"

    return app


if __name__ == "__main__":
    energies_df = pd.read_pickle("dash_data")
    app = get_app(energies_df)
    app.run_server(debug=True, reloader_type="stat")
