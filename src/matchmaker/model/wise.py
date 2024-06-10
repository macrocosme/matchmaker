# WISE Color-color cut for AGN from Jarret et al. (2011).
AGN_color_color_box = [(2.2, 1.7),                  (4.2, 1.7),
                                                    (4.2, 0.1 * 4.2 + 0.38),
                       (2.2, 0.1 * 2.2 + 0.38)]
def is_AGN(w1w2, w2w3):
    """
    Determines if an object is an Active Galactic Nucleus (AGN) based on the WISE color criteria.

    Args:
        w1w2 (float): The W1-W2 color index.
        w2w3 (float): The W2-W3 color index.

    Returns:
        bool: True if the object is classified as an AGN, False otherwise.
    """
    w1w2_cond = (w1w2 > (0.1 * w2w3 + 0.38)) and (w1w2 < 1.7)
    w2w3_cond = (w2w3 > 2.2) and (w2w3 < 4.2)
    return w1w2_cond & w2w3_cond
