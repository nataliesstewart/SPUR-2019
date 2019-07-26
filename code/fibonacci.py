def generate_basis_star(n,memo={}):
    if n == 2:
        return ["*p*"]
    elif n == 3:
        return ["*pp*"]
    
    basis = ["*"+"p"*(n-1)+"*"]
    for i in range(1,n-2):
        if n-i-1 not in memo.keys():
            memo[n-i-1] = generate_basis_star(n-i-1,memo)
        for elt in memo[n-i-1]:
            basis += ["*"+"p"*i+elt]
    return basis

def generate_basis_p(n):
    return [elt[:-1] for elt in generate_basis_star(n+1)]

def action(s,i,t):
    if s[:i+1] != t[:i+1] or s[i+2:] != t[i+2:]:
        return "0"

    if s[i:i+3] in ["*pp","pp*"]:
        return "a"

    if s[i:i+3] == "*p*":
        return "b"

    if s[i+1] != t[i+1]:
        return "d"

    if s[i+1] == "*":
        return "c"

    return "e"

def generate_matrices_star(n):
    """
    Generate the matrices for the fibonacci representation on n indices.
    """
    basis = generate_basis_star(n)
    out = [[[0 for i in range(len(basis))] for j in range(len(basis))] for k in range(n-1)]
    for i in range(n-1):
        for x,b in enumerate(basis):
            for y,c in enumerate(basis):
                out[i][x][y] = action(b,i,c)
    return out

def generate_matrices_p(n):
    """
    Generate the matrices for the fibonacci representation on n indices.
    """
    basis = generate_basis_p(n)
    out = [[[0 for i in range(len(basis))] for j in range(len(basis))] for k in range(n-1)]
    for i in range(n-1):
        for x,b in enumerate(basis):
            for y,c in enumerate(basis):
                out[i][x][y] = action(b,i,c)
    return out

for n in range(2,16):
    f = open("Fib m={} star".format(n),"w")
    f.write("\n/\n".join(["\n".join([",".join(row) for row in arr]) for arr in generate_matrices_star(n)]))
    f.close()
    
    f = open("Fib m={} p".format(n),"w")
    f.write("\n/\n".join(["\n".join([",".join(row) for row in arr]) for arr in generate_matrices_p(n)]))
    f.close()
    print("m={} completed.".format(n))
