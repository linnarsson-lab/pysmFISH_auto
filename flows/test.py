import numpy as np



def inner_function(idx, a):
    c = idx + a
    return c

def inner_function_2(idx, a, z):
    c = idx + a
    return c*1000

def runner(inner_function, d, **kwargs):
    c = inner_function(d[0], d[1])
    print(c)


if __name__ == "__main__":
    runner(inner_function_2, (20,10))