from egcd import egcd
def gcd(a : int, N : int) -> int:
    """
    returns the inverse element of a in modulo N
    """
    g, a_inv, _ = egcd(a,N)
    if g != 1:
        raise Exception('modular inverse does not exist') 
    else:
        if a_inv < 0:
            a_inv = a_inv % N
        return a_inv