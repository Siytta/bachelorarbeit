from egcd import egcd
def gcd(a, N):
    g, a_inv, _ = egcd(a,N)
    if g != 1:
        raise Exception('modular inverse does not exist') 
    else:
        if a_inv < 0:
            a_inv = a_inv % N
        return a_inv