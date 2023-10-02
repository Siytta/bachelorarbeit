from binToDez import binToDez, dezToBin, mod_exp
from gcd import gcd
from fractions import Fraction
from random import randint
from shor import Shor, Shor_sequential

import math

def kgV(a : int, b : int) -> int:
    """
    Returns the least common multiple of a and b
    """
    return abs(a * b) // math.gcd(a, b)

def fracPeriod(measuring: int,N: int, k: int) -> Fraction:
    """
    Returns the continued fractions of measuring/2**k limited to N
    """
    return Fraction(measuring/2**k).limit_denominator(N)

def isPeriod(a: int, N: int, denominator: int) -> bool:
    """
    Returns True if the denominator is the Periode of a in modulo N
    """
    return mod_exp(a, denominator, N) == 1

def passCriteria(a: int, N: int, fracMeasuring: Fraction) -> bool:
    """
    Checks the Criteria which is needed to extract the Factors of N from the result

    Parameters: 
    a:int
        The number which probably has a periodicity of fracMeasuring in modulo N
    N:int 
        The number of which the periodicity is tested 
    fracMeasuring:Fraction
        The number which dominator probably is the Periode of a in modulo N
    """
    return fracMeasuring.denominator%2==0 and mod_exp(a,fracMeasuring.denominator//2, N) != N-1

def calcFactors(a: int, N: int, fracMeasuring: Fraction) -> list[int,int]:
    """
    Calculates the factors of N from a and the periode which is in the denominator of fracMeasuring
    """
    return [math.gcd(a**(fracMeasuring.denominator//2)-1, N),math.gcd(a**(fracMeasuring.denominator//2)+1, N)]

def nonTrivialFactor(factors: list, N : int) -> bool:
    """
    Checks if atleast one of the two Elements of factors is a non trivial factor on N
    """
    return not ((factors[0] == 1 or factors[0] == N) and (factors[1] == 1 or factors[1] == N))

def isInMultiplePeriod(a : int ,N : int ,denominator: int, neightbore_range : int) -> bool:
    """
    Checks if close multiples of the denominator are the periode of a**p modulo N
    in the range of neightbore_range
    """
    for i in range(2,neightbore_range):
        if isPeriod(a, N, denominator*i):
            return True, Fraction(1, denominator*i)
    return False, 0

def isInLeastCommonMultiplePeriod(a: int,N: int, k: int , denominator:int , measurements_set: set) -> bool:
    """
    Checks if the least common multiple with one of the values of measurements_set 
    is the periode of a**p modulo N
    """
    for measurement in measurements_set:
        if isPeriod(a, N, kgV(denominator, measurement.frac.denominator)):
            return True, Fraction(1, kgV(denominator, measurement.frac.denominator))
    return False, 0

def hasNeightborePeriod(a:int, N:int, k:int, value:int, measurements_set:set, neightbore_range:int) -> bool: 
    """
    Checks if the close neightbore values of the messured value are the periode or a close multiple of the periode 
    """
    for i in range(1,neightbore_range):
        if value+i < 2**k:
            fracToTest = fracPeriod(value+i, N, k)
            if isPeriod(a , N, fracToTest.denominator):
                return True,  fracToTest
            condition, newFrac = isInMultiplePeriod(a, N, fracToTest.denominator, neightbore_range)
            if condition:
                return True, newFrac
            condition, newFrac = isInLeastCommonMultiplePeriod(a ,N ,k ,fracToTest.denominator ,measurements_set)
            if condition: 
                return True, newFrac
        if value-i > 0: 
            fracToTest = fracPeriod(value-i, N, k)
            if isPeriod(a , N, fracToTest.denominator):
                return True,  fracToTest
            condition, newFrac = isInMultiplePeriod(a, N, fracToTest.denominator, neightbore_range)
            if condition:
                return True, newFrac
            condition, newFrac = isInLeastCommonMultiplePeriod(a ,N ,k ,fracToTest.denominator ,measurements_set)
            if condition: 
                return True, newFrac
    return False, 0

def removeTrivial(factors : list, N : int) -> list: 
    """
    Given a list of two Elements one is atleast a nontrivial factor of N. 
    The function replaces the trivial factor of the elements in the lists with the other non-Trivial factor if not both elements are non-trivial factors 
    """
    if (factors[0] == 1 or factors[0] == N):
        factors[0] = N / factors[1]
    if (factors[1] == 1 or factors[1] == N):
        factors[1] = N / factors[0]
    return factors

class MeasuringValue:
    """
    Class to store the messured Value
    """
    def __init__(self, value, frac):
        self.value = value
        self.frac = frac
        
    def __str__(self):
        return f'Messurment: {self.value} with fraction: {self.frac}'

def Faktorisierungsalgorithmus(N: int, sequential: bool = True, m : int = -1, k : int = -1) -> list[int,int]:
    """
    The complete Algorithm to find the primefactors of a Number N.

    Parameters:
    N:int
        The number from which the prime factors are to be found
    sequential: bool:
        True if the Shor Algorithm with iterative Quanten-Phase-Estimation should be used.
        False if the Shor Algorithm with regular Quanten-Phase-Estimation should be used
    m: int
        Defines the approximationfactor of the Quantum-Fourier-Transformation
    k: int
        Defines the precision of the messurement(Bit count of the Result)
    """
    neightbore_range = 2**math.ceil(math.log2(N.bit_length()))
    if m == -1:
        m = N.bit_length() + 2
    if k == -1:
        k = 2*N.bit_length()
    while True:
        a = randint(2, N-1)
        if math.gcd(a, N) != 1:
            return [math.gcd(a, N), N//math.gcd(a, N)] #Zufällig Faktoren gefunden
        measured_fracs = set()
        quantum_circuit = 0
        while True:
            if sequential:
                measuring, quantum_circuit = Shor_sequential(a, N, quantum_circuit = quantum_circuit, m = m, k = k, number_shots = 1)
            else:
                measuring, quantum_circuit = Shor(a, N,quantum_circuit = quantum_circuit, m = m, k = k , number_shots = 1)
            measuring = MeasuringValue(int(list(measuring)[0], 2),fracPeriod(int(list(measuring)[0], 2),N, k))
            if not isPeriod(a, N, measuring.frac.denominator):
                condition, newFrac = isInMultiplePeriod(a, N, measuring.frac.denominator,neightbore_range)
                if condition:
                    measuring.frac = newFrac
                    break
                condition, newFrac = isInLeastCommonMultiplePeriod(a ,N ,k ,measuring.frac.denominator ,measurements_set)
                if condition:
                    measuring.frac = newFrac
                    break
                condition, newFrac = hasNeightborePeriod(a ,N , k,measuring.value , measurements_set, neightbore_range)
                if condition:
                    measuring.frac = newFrac
                    break
                measurements_set.add(measuring)
                continue
            break
        if not passCriteria(a, N, measuring.frac):
            continue
        factors = calcFactors(a, N, measuring.frac)
        if not nonTrivialFactor(factors, N):
            continue
        return removeTrivial(factors, N)