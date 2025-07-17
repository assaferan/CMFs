#!/usr/bin/env sage
"""
Sage script to compute q-expansions of Eisenstein cusp forms.

This script generates the q-expansion of the basis of the space of
Eisenstein series for given level, weight, and character.

To run, use the command:
    sage -python eisenstein_series_q_exp.py
"""

from sage.all import QQ, CyclotomicField, DirichletGroup, EisensteinForms
from sage.modular.dirichlet import DirichletCharacter
from lmfdb import db
from lmfdb.characters.TinyConrey import ConreyCharacter, get_sage_genvalues


# Additional parameters that might be useful
PRECISION = 10   # Number of terms in q-expansion
COND_BOUND = 19
WEIGHTS = [3, 4, 5]  # Weights to consider for Eisenstein series
OUTPUT_FILE = f"N1_{COND_BOUND}k{WEIGHTS[0]}-{WEIGHTS[-1]}lim{PRECISION}_trace0.sage.txt"

def num2letters(n):
    r"""
    Convert a number into a string of letters
    """
    if n <= 26:
        return chr(96+n)
    else:
        return num2letters(int((n-1)/26))+chr(97+(n-1) % 26)

def create_label(chi_orbit,weight, x):
    """
    Create a label for the Eisenstein series based on the character and Galois orbit code.

    Parameters:
    -----------
    chi : DirichletCharacter
        The Dirichlet character associated with the Eisenstein series.
    galois_orbit_code : str
        The Galois orbit code for the character.

    Returns:
    --------
    str
        A label string combining the character and Galois orbit code.
    """

    N,a = chi_orbit.split('.')
    N = int(N)  # Convert level N to integer
    k = weight
    label = f"{N}.{k}.{a}.E.{x}"
    return label

def trace_vector(eis_ser, prec=20):
    """
    Compute the trace vector of an Eisenstein series.

    Parameters:
    -----------
    eis_ser : EisensteinSeries
        The Eisenstein series for which to compute the trace vector

    Returns:
    --------
    list
        List of coefficients in the trace vector

    Notes:
    ------
    This function needs to be implemented using Sage's modular forms
    functionality.
    """
    return [a_i.trace() for a_i in eis_ser[0:prec+1]]


def get_coeffs_as_list_of_lists(eis_ser,K, prec=100):
    """
    Get the coefficients of the Eisenstein series as a list of lists.

    Parameters:
    -----------
    eis_ser : EisensteinSeries
        The Eisenstein series for which to get coefficients
    prec : int
        Number of terms in the q-expansion

    Returns:
    --------
    list
        List of coefficients in the q-expansion

    Notes:
    ------
    This function needs to be implemented using Sage's modular forms
    functionality.
    """
    if K.degree() == 1:
        coeffs = []
        for i in range(prec+1):
            coeffs.append([eis_ser[i]])

    coeffs = []
    for i in range(prec+1):
        coeff = K(eis_ser[i])
        coeff_as_vector = coeff.vector()
        coeffs.append(list(coeff_as_vector))
    return coeffs



def eisenstein_series_basis(level, weight, character=None, character_orbit=None, first=None, sage_zeta_order=None, precision=PRECISION, K=None):
    """
    Compute the q-expansion of the basis of Eisenstein series.

    Parameters:
    -----------
    level : int
        The level N of the modular forms
    weight : int
        The weight k of the modular forms
    character : DirichletCharacter or None
        The character (trivial if None)
    precision : int
        Number of terms in the q-expansion

    Returns:
    --------
    list
        List of q-expansions forming a basis of the Eisenstein subspace

    Notes:
    ------
    This function needs to be implemented using Sage's modular forms
    functionality. Key concepts to implement:
    - ModularForms(level, weight, character)
    - Eisenstein subspace
    - q-expansion computation
    """

    # TODO: Implement the actual computation
    # This will likely involve:
    # 1. Creating the modular forms space
    # 2. Getting the Eisenstein subspace
    # 3. Computing basis elements
    # 4. Getting q-expansions

    E = EisensteinForms(character,weight)
    E.set_precision(precision+1)
    new_eis_ser = E.new_eisenstein_series()
    trace_vecs = [trace_vector(eis_ser, precision) for eis_ser in new_eis_ser]
    sorted_trace_with_eis_ser = sorted(zip(trace_vecs, new_eis_ser), key=lambda x: x[0][1:])
    output = [(create_label(character_orbit, weight, num2letters(i+1)), trace_vec, get_coeffs_as_list_of_lists(eis_ser, K, precision)) for i, (trace_vec, eis_ser) in enumerate(sorted_trace_with_eis_ser)]
    output = [(character_orbit, first, weight, sage_zeta_order, label, trace_vec, eis_ser) for label, trace_vec, eis_ser in output]
    return output


def main():
    """
    Main function to compute and display the Eisenstein series basis.
    """
    print("Eisenstein Series Q-Expansion Generator")
    print("=" * 40)

    dirchar_table = db.char_dirichlet
    query = {
        'modulus': {'$gte' : 1,  '$lte': COND_BOUND}
        # 'is_primitive' : True,
            }
    payload = list(dirchar_table.search(query=query, projection=['conductor', 'modulus', 'first', 'label', 'order']))
    output = []
    for weight in WEIGHTS:
        for one_dir_char_orbit in payload:
            conductor = one_dir_char_orbit['conductor']
            modulus = one_dir_char_orbit['modulus']
            first = one_dir_char_orbit['first']
            dirchar_label = one_dir_char_orbit['label']
            order = one_dir_char_orbit['order']
            print(f"Conductor: {conductor}, Modulus: {modulus}, First: {first}, Label: {dirchar_label}")
            chi = ConreyCharacter(modulus, first)
            sage_zeta_order = chi.sage_zeta_order(order)
            genvalues_for_code = get_sage_genvalues(modulus, order, chi.genvalues, sage_zeta_order)
            K = CyclotomicField(sage_zeta_order)
            H = DirichletGroup(modulus, base_ring=K)
            M = H._module
            # if modulus == 9 and first == 1:
            #     import pdb; pdb.set_trace()  # Debugging point

            the_string = 'DirichletCharacter(H, M([{}]))'.format(
                    ','.join(str(val) for val in genvalues_for_code))
            sage_chi = eval(the_string)

            eis_ser_label_and_q_exp = eisenstein_series_basis(conductor, weight, sage_chi, dirchar_label, first, sage_zeta_order, PRECISION, K)
            output.append(eis_ser_label_and_q_exp)

    with open(OUTPUT_FILE, "w") as f:
        for eis_ser_one_dirchar_list in output:
            for character_orbit, first, weight, sage_zeta_order, label, trace_vec, eis_ser in eis_ser_one_dirchar_list:
                f.write(f"{character_orbit},{first},{weight},{PRECISION},{sage_zeta_order}," +
                        str(eis_ser).replace(' ', '') + "," +
                        str(trace_vec).replace(' ', '') + f",{label}\n")
    print(f"Results written to {OUTPUT_FILE}")


if __name__ == "__main__":
    # Run the main computation
    result = main()