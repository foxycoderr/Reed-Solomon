def i(n):
    return pow(n, -1, 251)

def m(a, b):
    return a*b%251

def plus(a, b):
    return (a+b)%251

#print(m(2, 245))
#print((7*250 + 8*2) % 251)

#print(pow(21, 2, 251))
print(plus(127, -233))