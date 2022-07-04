import astropy.units as u

def set_size(width=800, fraction=1, vertical=True):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Modified from https://jwalton.info/Embed-Publication-Matplotlib-Latex

    Parameters
    ----------
    width: float
            Width in pts
    fraction: float
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5 ** 0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio

    return (fig_width_in, fig_height_in) if not vertical else (fig_height_in, fig_width_in)

def set_fig_dims(direction, data_arr):
    if direction == 'horizontal':
        ncols = len(data_arr)
        nrows = 1
    elif direction == 'vertical':
        ncols = 1
        nrows = len(data_arr)

    return ncols, nrows
