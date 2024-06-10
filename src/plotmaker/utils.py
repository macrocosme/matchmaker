import astropy.units as u

GOLDEN_RATIO = (5 ** 0.5 - 1) / 2
# WIDTH=3.54

WIDTH_twocols = (7.1 - 0.3125)/2  # From A&A latex template
WIDTH = WIDTH_twocols  # Set default fot the two column A&A format
WIDTH_onecol = 7.1
DPI=600

def set_size(width=3.54, dpi=150, fraction=1, vertical=True):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Modified from https://jwalton.info/Embed-Publication-Matplotlib-Latex

    Parameters
    ----------
    width: astropy.units.cm
            Width in astropy.units.cm
    dpi: int
            dot per inch
    fraction: float
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Figure in inches
    fig_width_in = width * fraction * dpi
    fig_height_in = fig_width_in * GOLDEN_RATIO

    return (fig_width_in, fig_height_in) if not vertical else (fig_height_in, fig_width_in)

def set_fig_dims(direction, data_arr):
    if direction == 'horizontal':
        ncols = len(data_arr)
        nrows = 1
    elif direction == 'vertical':
        ncols = 1
        nrows = len(data_arr)

    return ncols, nrows
