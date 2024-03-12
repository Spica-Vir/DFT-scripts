#!/usr/bin/env python

###############################################################################
# To read the output eigenvalues of the overlap matrix from CRYSTAL .out file.#
# Output: The minimum eigenvalue and its corresponding k point; Whether there #
# is negative eigenvalue; Whether the basis is linear dependent.              #
# ----------------------------                                                #
# by Spica. Vir., ICL, Nov. 20 - 21                                           #
# spica.h.zhou@gmail.com                                                      #
###############################################################################

def readeigs(eigs):
    countline = 0
    label = 'S(K) EIGENV - K'
    mi_eigs = []
    mx_eigs = []
    while countline < len(eigs) - 1:
        title = ' '.join(eigs[countline].strip().split()[0:4])
        if title == label:
            countline += 1
            eigs_k = []
            while len(eigs[countline].strip().split()) and \
                    eigs[countline].strip().split()[0][0] != 'T':
                eigs_k.extend(eigs[countline].strip().split())
                countline += 1

            eigs_k = list(map(float, eigs_k))
            mi_eigs.append(min(eigs_k))
            mx_eigs.append(max(eigs_k))
        else:
            countline += 1

    mi_eig = min(mi_eigs)
    mi_k = mi_eigs.index(min(mi_eigs))
    mx_eig = max(mx_eigs)
    mx_k = mx_eigs.index(max(mx_eigs))
    return mi_eig, mi_k, mx_eig, mx_k


filename = input('The name of .out file = ')
file = open(filename, "r")
eigs = file.readlines()
file.close()

mi_eig, mi_k, mx_eig, mx_k = readeigs(eigs)

if mi_k < 0:
    print('Warning: Basis set linear dependency detected.')

print('The minimum eigenvalue = %8.6f at K = %d' % (mi_eig, mi_k))
print('The maximum eigenvalue = %8.6f at K = %d' % (mx_eig, mx_k))
