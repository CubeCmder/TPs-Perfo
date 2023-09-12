
def knots2fps(v, reverse=False):
    """
    Convert knots to feet per second. If reverse is set to True, does the opposite operation.

    :param v: Value to convert
    :param reverse: True if the opposite operation is required, i.e. [f/s] to [knots]

    :return: Converted value.
    """

    if reverse is False:
        return v*1.68781
    else:
        return v/1.68781