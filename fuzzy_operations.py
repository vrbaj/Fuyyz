"""
.. versionadded:: 0.2

Module for operations with fuzzy sets. Namely complement, intersection (T-Norms) and union (S-Norms)

Contents of this page:

.. contents::
    :local:
    :depth: 1


Fuzzy complement, union (S-norms) and intersection (T-norms)
============================================================


"""


def complement_standard(membership_degree):
    complement = 1.0 - membership_degree
    return complement


def complement_sugeno(membership_degree, lambda_value):
    if lambda_value <= -1:
        print("sugeno_complement:parameter lambda out of range, it should be (-1,inf)")
    complement = (1 - membership_degree) / (1 + lambda_value * membership_degree)
    return complement


def complement_yager(membership_degree, w_value):
    """
    Yager
    """
    if w_value <= 0:
        print("yager_complement: parameter w out of range, it should be (0, inf)")
    complement = (1 - membership_degree ** w_value) ** (1 / w_value)
    return complement


def s_norm_max(membership_degreee_a, membership_degree_b):
    s_norm_value = max(membership_degreee_a, membership_degree_b)
    return s_norm_value


def s_norm_dombi(membership_degree_a, membership_degree_b, lambda_value):
    if lambda_value <= 0:
        print("parameter lambda out of bounds, it should be (0, inf)")
    s_norm_value = 1 / (1 + ((1 / membership_degree_a - 1) ** (-1 * lambda_value) +
                        (1 / membership_degree_b - 1) ** (-1 * lambda_value)) ** (-1 / lambda_value))
    return s_norm_value


def s_norm_dubois_prade(membership_degree_a, membership_degree_b, alpha):
    if not(0 <= alpha <= 1):
        print("parameter alpha out of bounds, it should be [0,1]")
    s_norm_value = (membership_degree_a + membership_degree_b - membership_degree_a * membership_degree_b -
                    min(membership_degree_a, membership_degree_b, 1 - alpha)) / \
                   (max(1 - membership_degree_a, 1 - membership_degree_b, alpha))
    return s_norm_value


def s_norm_yager(membership_degree_a, membership_degree_b, w):
    """
    :param float membership_degree_a: membership degree to set A
    :param float membership_degree_b: membership degree to set B
    :param float w: w parameter
    :return: membership degree of union
    """
    if w <= 0:
        print("parameter w is out of bounds, it should be (0, inf)")
    s_norm_value = min(1.0, ((membership_degree_a ** w) + (membership_degree_b ** w)) ** (1 / w))
    return s_norm_value


def s_norm_drastic(membership_degree_a, membership_degree_b):
    s_norm_value = 1
    if membership_degree_a == 0:
        s_norm_value = membership_degree_b
    if membership_degree_b == 0:
        s_norm_value = membership_degree_a
    return s_norm_value


def s_norm_einstein(membership_degree_a, membership_degree_b):
    s_norm_value = (membership_degree_a + membership_degree_b) / (1 + membership_degree_a * membership_degree_b)
    return s_norm_value


def s_norm_algebraic(membership_degree_a, membership_degree_b):
    s_norm_value = membership_degree_a + membership_degree_b - membership_degree_a * membership_degree_b
    return s_norm_value


def t_norm_min(membership_degree_a, membership_degree_b):
    s_norm_value = min(membership_degree_a, membership_degree_b)
    return s_norm_value


def t_norm_dombi(membership_degree_a, membership_degree_b, lambda_value):
    if lambda_value <= 0:
        print("parameter lambda out of bounds, it should be (0, inf)")
    t_norm_value = 1 / (1 + ((1 / membership_degree_a - 1) ** lambda_value +
                        (1 / membership_degree_b - 1) ** lambda_value) ** (1 / lambda_value))
    return t_norm_value


def t_norm_dubois_prade(membership_degree_a, membership_degree_b, alpha):
    if not(0 <= alpha <= 1):
        print("parameter alpha out of bounds, it should be [0, 1]")
    t_norm_value = membership_degree_a * membership_degree_b / max(membership_degree_a, membership_degree_b, alpha)
    return t_norm_value


def t_norm_yager(membership_degree_a, membership_degree_b, w):
    if w <= 0:
        print("parameter w is out of bounds, it should be (0, inf")
    t_norm_value = 1 - min(1.0, ((1 - membership_degree_a) ** w + (1 - membership_degree_b) ** w) ** (1 / w))
    return t_norm_value


def t_norm_drastic(membership_degree_a, membership_degree_b):
    t_norm_value = 0
    if membership_degree_a == 1:
        t_norm_value = membership_degree_b
    if membership_degree_b == 1:
        t_norm_value = membership_degree_a
    return t_norm_value


def t_norm_einstein(membership_degree_a, membership_degree_b):
    t_norm_value = membership_degree_a * membership_degree_b / (2 - (membership_degree_a + membership_degree_b -
                                                                     membership_degree_a * membership_degree_b))
    return t_norm_value


def t_norm_product(membership_degree_a, membership_degree_b):
    t_norm_value = membership_degree_a * membership_degree_b
    return t_norm_value


# testing_membership_degree = 0.9
# complement_value = complement_yager(testing_membership_degree, 0.3)
# print(complement_value)

testing_membership_degree_a = 0.7
testing_membership_degree_b = 1
complement_value = t_norm_product(testing_membership_degree_a, testing_membership_degree_b)
print(complement_value)
