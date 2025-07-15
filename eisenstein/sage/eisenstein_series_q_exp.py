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


WEIGHT = 3       # Weight k

# Additional parameters that might be useful
PRECISION = 20   # Number of terms in q-expansion

OUTPUT_FILE = f"eisenstein_series_qexp_weight_{WEIGHT}.txt"

def num2letters(n):
    r"""
    Convert a number into a string of letters
    """
    if n <= 26:
        return chr(96+n)
    else:
        return num2letters(int((n-1)/26))+chr(97+(n-1) % 26)

def create_label(chi_orbit, x):
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
    k = WEIGHT
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
    return [a_i.trace() for a_i in eis_ser[1:prec]]


def eisenstein_series_basis(level, weight, character=None, character_orbit=None, precision=10):
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
    E.set_precision(precision)
    new_eis_ser = E.new_eisenstein_series()
    trace_vecs = [trace_vector(eis_ser, precision) for eis_ser in new_eis_ser]
    sorted_trace_with_eis_ser = sorted(zip(trace_vecs, new_eis_ser), key=lambda x: x[0])
    output = [(create_label(character_orbit, num2letters(i+1)), eis_ser) for i, (_, eis_ser) in enumerate(sorted_trace_with_eis_ser)]
    output = {label: eis_ser for label, eis_ser in output}
    return output


def main():
    """
    Main function to compute and display the Eisenstein series basis.
    """
    print("Eisenstein Series Q-Expansion Generator")
    print("=" * 40)

    dirchar_table = db.char_dirichlet
    query = {
        'modulus': {'$gte' : 1,  '$lte': 20},
        'is_primitive' : True,
            }
    payload = dirchar_table.search(query=query, projection=['conductor', 'first', 'label', 'order'])
    output = []
    for one_dir_char_orbit in payload:
        conductor = one_dir_char_orbit['conductor']
        first = one_dir_char_orbit['first']
        label = one_dir_char_orbit['label']
        order = one_dir_char_orbit['order']
        print(f"Conductor: {conductor}, First: {first}, Label: {label}")
        chi = ConreyCharacter(conductor, first)
        sage_zeta_order = chi.sage_zeta_order(order)
        genvalues_for_code = get_sage_genvalues(conductor, order, chi.genvalues, sage_zeta_order)
        H = DirichletGroup(conductor, base_ring=CyclotomicField(sage_zeta_order))
        M = H._module

        the_string = 'DirichletCharacter(H, M([{}]))'.format(
                ','.join(str(val) for val in genvalues_for_code))
        sage_chi = eval(the_string)


        eis_ser_label_and_q_exp = eisenstein_series_basis(conductor, WEIGHT, sage_chi, label, PRECISION)
        output.append(eis_ser_label_and_q_exp)

    with open(OUTPUT_FILE, "w") as f:
        for eis_ser_dict in output:
            for label, eis_ser in eis_ser_dict.items():
                f.write(f"{label},{eis_ser}\n")
    print(f"Results written to {OUTPUT_FILE}")


if __name__ == "__main__":
    # Run the main computation
    result = main()