def binToDez(bin: list[int]) -> int:
    dez = 0
    index = 0
    for e in range(len(bin)):
        dez += (2**e)* bin[index]
        index+=1
    return dez

def dezToBin(dez: int, length: int = -1) -> list[int]:
    bin = list(reversed([int(i) for i in list('{0:0b}'.format(dez))]))
    if length != -1:
        while len(bin) < length:
            bin = bin + [0]
    return bin

    
def mod_exp(base : int , exponent : int, modulus : int) -> int:
    """
    returns base**exponent % modulus 
    """
    result = 1
    base = base % modulus
    while exponent > 0:
        if exponent % 2 == 1:
            result = (result * base) % modulus
        exponent = exponent >> 1
        base = (base * base) % modulus
    return result   

