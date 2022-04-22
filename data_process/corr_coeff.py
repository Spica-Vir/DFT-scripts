#!/usr/bin/env python
# coding: utf-8

# A corr_coeff function to calculate the pairwise correlation coefficient
# Requirement: pandas package, imported as pd
# By Spica. Vir., ICL, Apr. 22, 22. spica.h.zhou@gmail.com

import pandas as pd


def corr_coeff(varA, varB, varA_name='varA', varB_name='varB',
               method='spearman'):
    """
    function corr_coeff, to calculate the correlation coefficient between 2
    variables.
    varA: list, the list of values of variable A
    varB: list, the list of values of variable B
    varA_name: str, the name of variable A, used for indexing
    varB_name: str, the name of variable B, used for indexing
    method: str, correlation methods: 'pearson', 'kendall', 'spearman'
    """
    if len(varA) != len(varB):
        print('Error: variable vectors should have the same dimension.')
        return

    pairs = [(varA[i], varB[i], ) for i in range(len(varA))]
    correlated_data = pd.DataFrame(pairs, columns=[varA_name, varB_name])
    corr_coeff = correlated_data.corr(method=method)

    return corr_coeff
