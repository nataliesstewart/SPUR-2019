import cmath
import numpy
import crossingless
import time
import csv
import math

def make_matrix(n,q,restriction = 1):
    """

    Create the matrix for the direct sum (1 + t_1) oplus ... oplus (1 + t_{i - 1 - restriction}) for parameter q
    This has trivial kernel iff the restriction to S_{n - restriction} is irreducible.
    
    KEYWORD ARGUMENTS:
    n -- parameter from S_n
    q -- parameter of the hecke algebra (T_i - q)(T_i + 1) = 0
    restriction -- the number of simple transpositions ignored.
    """
    assert n % 2 == 0

    # generate a basis; basisraw has basis_element elements, basis has module_element
    basisraw = crossingless.generate_basis(n)
    basis = [crossingless.module_element(n,{b:1},q) for b in basisraw]

    C = len(basis)
    m = numpy.zeros(C)
    # i iterates over the simple transpositions considered
    for i in range(n - restriction - 1):
        # generate the matrix by applying the transposition to each basis vector and taking the resulting coefficients
        outs = [b.mult_generator(i+1) for b in basis]
        outsmat = numpy.matrix([[out.coeff(b) for out in outs] for b in basisraw])
        m = numpy.vstack([m,outsmat])
    
    # prune zero rows from matrix
    m = m[[i for i, x in enumerate(m) if x.any()]]
    return m

def representation(n,q):
    """
    Create the matrices spelling out the representation
    """
    # generate a basis; basisraw has basis_element elements, basis has module_element
    basisraw = crossingless.generate_basis(n)
    basis = [crossingless.module_element(n,{b:1},q) for b in basisraw]

    C = len(basis)
    m = []
    # i iterates over the simple transpositions considered
    for i in range(n - 1):
        # generate the matrix by applying the transposition to each basis vector and taking the resulting coefficients
        outs = [b.mult_generator(i+1) for b in basis]
        outsmat = numpy.matrix([[out.coeff(b) for out in outs] for b in basisraw])
        m = m + [outsmat]
    
    return m
    

def kernel_dimension(n,q,restriction=1):
    """
    Find the dimension of the kernel of the matrix found above.
    """
    m = make_matrix(n,q,restriction)
    return m.shape[1] - numpy.linalg.matrix_rank(m)

def test_modular_kernels(nmax,nmin=4,restrictions=[1],debug=False,write=None):
    """
    Test all modular crossingless matchings representations (e >= 3) for nmin <= n <= nmax and the corresponding restriction.
    
    KEYWORD ARGUMENTS:
    nmax -- highest n tested
    nmin -- lowest n tested
    restrictions -- tests at S_{n-r} for each r.
    debug -- boolean flag enables or disables terminal output with times/values
    write -- csvwriter object to write results to
    """
    assert nmax % 2 == 0 and nmin % 2 == 0 and nmax >= nmin
    
    allkernels = []
    for restriction in restrictions:
        # write row of e's at the top of the file
        if write is not None:
            write.writerow(["Kernel ranks of restriction to S_(n - {})".format(restriction)])
            write.writerow(["n \\ e"] + list(range(3,nmax + 1)))
        
        # non-tested e's are blank
        kernels = [["" for j in range(nmax - 2)] for i in range((nmax - nmin)//2 + 1)]        

        if debug: # times for trials are just clock-time differences
            oldtime = time.process_time()
        
        # ignore n which are too small for the restriction
        restriction_rounded = math.ceil((restriction+2)/2.)*2 
        nmin_corrected = nmin if (restriction <= nmin - 2) else restriction_rounded
        
        for n in range(nmin_corrected,nmax + 1,2): # n ranges over the evens in [nmin,nmax]
            for e in range(3,n+1): # e ranges over the integers in [e,n]
                # make the primitive kth root exp(2 pi i/k) for k = e (confusing notation, I know)
                q = cmath.exp(2*cmath.pi*cmath.sqrt(-1)/e)
             
                dim = kernel_dimension(n,q,restriction)
                kernels[(n - nmin_corrected)//2][e - 3] = dim
                
                if debug:
                    newtime = time.process_time()
                    print("Case n = {:2d} \t e = {:2d} \t r = {:2d} has kernel dimension {} \t and took {:4f} seconds.".format(n,e,restriction,dim,newtime - oldtime))
                    oldtime = newtime
            
            if write is not None:
                write.writerow([n] + kernels[(n-nmin_corrected)//2])
        allkernels = allkernels + [kernels]
    return allkernels

# what follows is an example output, writing to file "Heuristic_data.csv" the kernels of
# S_{n - r} for r from 0 to high enough for nontrivial kernels.
"""
nmin = 4
nmax = 12
debug = True
restrictions = [0,1,2,3,4,5,6]
write = csv.writer(open("Heuristic_data.csv","w"))
test_modular_kernels(nmax,nmin,restrictions,debug,write)
"""
"""
# another test simply printing kernels for a particular q and restriction
nmin = 4
nmax = 16
q = 1
restriction = 1
oldtime = time.process_time()
for n in range(nmin,nmax+1,2):
    dim = kernel_dimension(n,q,restriction)
    newtime = time.process_time()
    print("Semisimple case n = {:2d} \t has kernel dmension {} \t and took {:4f} seconds.".format(n,dim,newtime - oldtime))
    oldtime = newtime
"""

q = 1
for n in range(4,17,2):
    write = csv.writer(open("Representation n={}.csv".format(n),"w"))
    rep = representation(n,q)
    rep = [arr.astype(int) for arr in rep]
    for arr in rep:
        write.writerows(arr)
        write.writerow(["/"])
    print("n={} completed.".format(n))
