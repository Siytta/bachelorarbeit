def binToDez(bin: list[int]):
    dez = 0
    index = 0
    for e in reversed(list(range(len(bin)))):
        dez += (2**e)* bin[index]
        index+=1
    return dez

