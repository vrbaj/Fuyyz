import numpy as np


def get_membership_degree(point, fuzzy_set_type, fuzzy_set_params):
    membership_degree = 0
    if fuzzy_set_type == "TRIANGULAR":
        membership_degree = get_triangular_mf_degree(point, fuzzy_set_params)
    elif fuzzy_set_type == "GAUSSIAN":
        membership_degree = get_gaussian_mf_degree(point, fuzzy_set_params)
    elif fuzzy_set_type == "GAUSSIAN2":
        membership_degree = get_gaussian2_mf_degree(point, fuzzy_set_params)
    elif fuzzy_set_type == "TRAPEZOIDAL":
        membership_degree = get_trapezoidal_mf_degree(point, fuzzy_set_params)
    elif fuzzy_set_type == "R":
        membership_degree = get_r_mf_degree(point, fuzzy_set_params)
    elif fuzzy_set_type == "L":
        membership_degree = get_l_mf_degree(point, fuzzy_set_params)
    elif fuzzy_set_type == "BELL":
        membership_degree = get_bell_mf_degree(point, fuzzy_set_params)
    elif fuzzy_set_type == "SIGMOID":
        membership_degree = get_sigmoid_mf_degree(point, fuzzy_set_params)
    return membership_degree


def get_triangular_mf_degree(point, params):
    """
    This function returns the membership degree of point to fuzzy set specified by triangular  member ship function.
     Triangular membership function is specified by parameters in params

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of triangular membership function [lower_limit, center, upper_limit] or
                 in symmetrical case [distance_from_center_to_bounds, center]

    Returns:
        float -- value of membership degree to triangular function

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
        elif center < point < upper_limit:
            membership_degree = (upper_limit - point)/(upper_limit - center)
        else:
            membership_degree = 0

    return membership_degree


def get_gaussian_mf_degree(point, params):
    """
    This function returns the membership degree of point to fuzzy set specified by gaussian  member ship function.
     Gaussian membership function is specified by parameters in params. Gaussian function is implemented as
     math.exp((-0.5*(point - center)**2)/(sigma**2))

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of gaussian membership function [center, sigma]

    Returns:
        float -- value of membership degree to gaussian function

    """
    membership_degree = 0
    if params.size <= 1 or params.size > 2:
        print('get_gaussian_mf_value: invalid number of parameters')
    else:
        center = params[0]
        sigma = params[1]
        membership_degree = np.exp((-0.5*(point - center)**2)/(sigma**2))
    return membership_degree


def get_gaussian2_mf_degree(point, params):
    """
    This function returns the membership degree of point to fuzzy set specified by 2 gaussian  member ship function.
     Gaussian membership functions are specified by parameters in params. The value of membership degree between two
     centers of gaussian functions is always 1. Gaussian function is implemented as
     np.exp((-0.5*(point - center)**2)/(sigma**2))

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of gaussian2 membership function [center1, sigma1, center2, sigma2]

    Returns:
        float -- value of membership degree to gaussian2 function

    """
    membership_degree = 0
    if params.size <= 3 or params.size > 4:
        print('get_gaussian2_mf_value: invalid number of parameters')
    else:
        center1 = params[0]
        sigma1 = params[1]
        center2 = params[2]
        sigma2 = params[3]
        if point < center1:
            membership_degree = np.exp((-0.5*(point - center1)**2)/(sigma1**2))
        elif center1 <= point <= center2:
            membership_degree = 1
        elif point > center2:
            membership_degree = np.exp((-0.5*(point - center2)**2)/(sigma2**2))
    return membership_degree


def get_trapezoidal_mf_degree(point, params):
    """
    This function returns the membership degree of point to fuzzy set specified by trapezoidal  member ship function.
     Trapezoidal membership function is specified by parameters in params.

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of trapezoidal membership function [a, b, c, d]

    Returns:
        float -- value of membership degree to trapezoidal function

    """
    membership_degree = 0
    if params.size <= 3 or params.size > 4:
        print("get_trapezoidal_mf_degree: invalid number of parameters")
    else:
        a = params[0]
        b = params[1]
        c = params[2]
        d = params[3]
        if point < a or point > d:
            membership_degree = 0
        elif a <= point <= b:
            membership_degree = (x - a) / (b - a)
        elif b <= point <= c:
            membership_degree = 1
        elif c <= point <= d:
            membership_degree = (d - point) / (d - c)
    return membership_degree


def get_r_mf_degree(point, params):
    """
    This function returns the membership degree of point to fuzzy set specified by R  member ship function.
     R membership function is specified by parameters in params.

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of R membership function [c, d]

    Returns:
        float -- value of membership degree to R function

    """
    membership_degree = 0
    if params.size <= 1 or params.size > 2:
        print("get_r_mf_degree: invalid number of parameters")
    else:
        c = params[0]
        d = params[1]
        if point > d:
            membership_degree = 0
        elif c <= point <= d:
            membership_degree = (d - point) / (d - c)
        elif point < c:
            membership_degree = 1
    return membership_degree


def get_l_mf_degree(point, params):
    """
    This function returns the membership degree of point to fuzzy set specified by L  member ship function.
     L membership function is specified by parameters in params.

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of L membership function [a, b]

    Returns:
        float -- value of membership degree to L function

    """
    membership_degree = 0
    if params.size <= 1 or params.size > 2:
        print("get_l_mf_degree: invalid number of parameters")
    else:

        a = params[0]
        b = params[1]
        if point > b:
            membership_degree = 1
        elif a <= point <= b:
            membership_degree = (point - a) / (b - a)
        elif point < a:
            membership_degree = 0
    return membership_degree


def get_bell_mf_degree(point, params):
    """
    This function returns the membership degree of point to fuzzy set specified by bell-shaped  member ship function.
     Bell-shaped membership function is specified by parameters in params. Bell-shaped function is implemented as
     1/(1+abs((point-center)/a)**2b

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of bell-shaped membership function [a, b, center]

    Returns:
        float -- value of membership degree to L function

    """
    membership_degree = 0
    if params.size <= 2 or params.size > 3:
        print("get_bell_mf_degree: invalid number of parameters")
    else:

        a = params[0]
        b = params[1]
        c = params[2]
        membership_degree = 1/(1+(abs((point - c) / a)**(2*b)))
    return membership_degree


def get_sigmoid_mf_degree(point, params):
    """
    This function returns the membership degree of point to fuzzy set specified by sigmoidal  member ship function.
     Sigmoidal membership function is specified by parameters in params. This function is implemented as
     1/(1+exp(-a*(point-center)).

    Args:
        point (number):  point in which we want to estimate membership degree
        params (numpy array):  specification of bell-shaped membership function [center, a]

    Returns:
        float -- value of membership degree to L function

    """
    membership_degree = 0
    if params.size <= 1 or params.size > 2:
        print("get_sigmoid_mf_degree: invalid number of parameters")
    else:

        center = params[0]
        a = params[1]
        membership_degree = 1/(1+np.exp(-a*(point - center)))
    return membership_degree
x = 2.01

test_fuzzy_set_params = np.array([2, 4])
# test_fuzzy_set_type = "TRIANGULAR"
# test_fuzzy_set_type = "GAUSSIAN"
# test_fuzzy_set_type = "TRAPEZOIDAL"
# test_fuzzy_set_type = "R"
# test_fuzzy_set_type = "L"
# test_fuzzy_set_type = "GAUSSIAN2"
test_fuzzy_set_type = "SIGMOID"
mf = get_membership_degree(x, test_fuzzy_set_type, test_fuzzy_set_params)
print(mf)
