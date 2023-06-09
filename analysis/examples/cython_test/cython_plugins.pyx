import sys


# A
cdef double a_adder = -1

def a_wrapper():
    print('I know about "a"')
    return 'a', a_setup, a_call

cdef a_setup(double adder):
     global a_adder
     print('Setting up "a"')
     a_adder = adder

cdef double a_call(double i):
    print('a:', i)
    return i + a_adder


# B
cdef b_multiplier = 3

def b_wrapper():
    return 'b', b_setup, b_call

def b_setup(multiplier):
    global b_multiplier
    b_multiplier = multiplier

def b_call(i):
    print('b:', i)
    return i * b_multiplier


def c(i):
    print('c:', i)
    return i / 3


def run():
    wrappers = [w for w in dir(sys.modules[__name__]) if w.endswith('_wrapper')]
    names, setup_funcs, call_funcs = zip(*[getattr(sys.modules[__name__], f)() for f in wrappers])
    print(names)

    for f_name, f_setup in zip(names, setup_funcs):
        print('Setting up ',f_name)
        f_setup(2)

    for n in range(20):
        for f_name, f_call in zip(names, call_funcs):
            print('Calling', f_name)
            print(f_call(n))
