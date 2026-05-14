from matplotlib import colormaps
import numpy as np
import plotly.graph_objects as go


def _to_plotly_colorscale(colormap: str, steps: int = 256):
    cmap = colormaps[colormap]
    return [
        [
            index / (steps - 1),
            f"rgb({int(red * 255)}, {int(green * 255)}, {int(blue * 255)})",
        ]
        for index, (red, green, blue, _) in enumerate(cmap(np.linspace(0, 1, steps)))
    ]


def _normalize_pae(pae):
    pae_array = np.asarray(pae, dtype=float)
    if pae_array.ndim != 2:
        raise ValueError("pae must be a 2D array-like object")
    return pae_array


def plot_pae(
    pae: np.ndarray | list,
    width=1000,
    height=1000,
    colormap="coolwarm",
    vmin=0,
    vmax=31.75,
):
    pae = _normalize_pae(pae)
    fig = go.Figure(
        data=go.Heatmap(
            z=pae,
            colorscale=_to_plotly_colorscale(colormap),
            zmin=vmin,
            zmax=vmax,
            colorbar={"title": "PAE (Å)"},
        )
    )
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title="Estimated Residue index",
        yaxis_title="Aligned Residue index",
    )
    fig.update_yaxes(autorange="reversed")
    return fig
