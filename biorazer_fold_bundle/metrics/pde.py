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


def _normalize_pde(pde):
    pde_array = np.asarray(pde, dtype=float)
    if pde_array.ndim != 2:
        raise ValueError("pde must be a 2D array-like object")
    return pde_array


def plot_pde(
    pde: np.ndarray | list,
    width=1000,
    height=1000,
    colormap="coolwarm",
    vmin=0,
    vmax=20,
):
    pde = _normalize_pde(pde)
    fig = go.Figure(
        data=go.Heatmap(
            z=pde,
            colorscale=_to_plotly_colorscale(colormap),
            zmin=vmin,
            zmax=vmax,
            colorbar={"title": "PDE (Å)"},
        )
    )
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title="Residue index",
        yaxis_title="Residue index",
    )
    fig.update_yaxes(autorange="reversed")
    return fig
