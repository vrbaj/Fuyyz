"""
.. versionadded:: 0.2

Module for operations with fuzzy sets. Namely complement, intersection (T-Norms) and union (S-Norms)

Contents of this page:

.. contents::
    :local:
    :depth: 1


Fuzzy complement
=============================
"""

import numpy as np


def standard_complement(membership_degree):
    complement = 1.0 - membership_degree
    return complement

def sugeno_complement(membership_degree, lambda_value):
    if lambda_value <= -1:
        print("sugeno_complement:parameter lambda out of range, it should be (-1,inf)")
    complement = (1 - membership_degree) / (1 + lambda_value * membership_degree)
    return complement


def yager_complement(membership_degree, w_value):
    if w_value <= 0:
        print("yager_complement: parameter w out of range, it should be (0, inf)")
    complement = (1 - membership_degree ** w_value) ** (1 / w_value)
    return complement


"""
Fuzzy union
=============================
"""

testing_membership_degree = 0.9
complement_value = yager_complement(testing_membership_degree,0.3)
print(complement_value)