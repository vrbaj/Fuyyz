import numpy as np


def get_membership_degree(point, fuzzy_set_type, fuzzy_set_params):
    membership_degree = 0
    if fuzzy_set_type == "TRIANGULAR":
        membership_degree = get_triangular_mf_degree(point, fuzzy_set_params)
    return membership_degree


def get_triangular_mf_degree(point, params):
    """
    This function return the membership degree of point to fuzzy set specified by triangular  member ship function.
     Triangular membership function is specified by parameters in params

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of triangular membership function [lower_limit, center, upper_limit] or
                 in symmetrical case [distance_from_center_to_bounds, center]

    Returns:
        float -- value of membership to triangular function

    """
    membership_degree = 0
    if params.size <= 1 or params.size > 3:
        print('get_triangular_mf_value: invalid number of parameters')
    else:
        center = params[1]
        if params.size == 2:
            lower_limit = params[1] - params[0]
            upper_limit = params[1] + params[0]
        elif params.size == 3:
            lower_limit = params[0]
            upper_limit = params[2]
        if lower_limit < point <= center:
            membership_degree = (point - lower_limit)/(center - lower_limit)
            print('jsem tu')
        elif center < point < upper_limit:
            membership_degree = (upper_limit - point)/(upper_limit - center)
            print('jsem tu podruhe')
        else:
            membership_degree = 0

    return membership_degree

x = -1
fuzzy_set_params = np.array([0, 1, 2])
fuzzy_set_type = "TRIANGULAR"
mf = get_membership_degree(x, fuzzy_set_type, fuzzy_set_params)
print(mf)
