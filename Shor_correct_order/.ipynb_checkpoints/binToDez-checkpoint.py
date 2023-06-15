def binToDez(bin: list[int]) -> int:
    dez = 0
    index = 0
    for e in reversed(list(range(len(bin)))):
        dez += (2**e)* bin[index]
        index+=1
    return dez

def dezToBin(dez: int, length: int = -1) -> list[int]:
    bin = [int(i) for i in list('{0:0b}'.format(dez))]
    if length != -1:
        while len(bin) < length:
            bin = [0] + bin
    return bin
    

