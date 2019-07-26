from cmath import sqrt

# n: number of nodes (note that this must be even)
# match: dictionary sending node to matching node
# for example, {1:4,4:1,2:3,3:2} is the S4 rainbow
#
#


class basis_element:
    def __init__(self, n, match, r=0):
        assert n % 2 == 0, 'n must be even'
        
        self.n = n
        self.r = r
        assert len(match.keys()) == n, 'must have n nodes, do not specify r'
        for i in match.keys():
            assert match[i] != i, 'an element cannot be paired with itself'
            assert i == match[match[i]], 'matchings are not symmetric'
            for j in range(min(i, match[i])+1, max(i, match[i])):
                assert min(i, match[i]) < match[j] and max(
                    i, match[i]) > match[j], 'matchings cross'
        
        # anchors match to -1
        # note: this does not check validity.        
        for i in range(1,n + r + 1): 
            if i not in match.keys():
                match[i] = i
        self.match = match

    # returns dictionary of matches nodes
    def get_matches(self):
        dictt = {}
        for i in self.match.keys():
            dictt[i] = self.match[i]
        return dictt

    # returns matching node to i
    def pair(self, i):
        return self.match[i]

    # creates tex to call the Matchings macro on the given basis elt
    def output_latex(self):
        prematching = []
        for key in self.match.keys():
            if key < self.match[key]:
                prematching += ["{}/{}".format(key, self.match[key])]

        matching = "\\Matching{{{}}}{{{}}}".format(
            self.n, ",".join(prematching))
        return matching

    def __eq__(self, i):
        return i.get_matches() == self.get_matches()

    # so that you can use the basis elements to index in a dictionary
    def __hash__(self):
        a = 0
        for i in range(self.n):
            a = a+self.n**i*self.match[i+1]
        return a

    def __str__(self):
        return str(self.match)


# vecpairs: dictionary sending basis element to coefficient
# for example, the element (rainbow+two humps) in S4 would be implemented
# using a dictionary sending the rainbow basis element to 1, and the two humps
# element to 1
class module_element:
    def __init__(self, n, vecpairs, q=1, r=0):
        assert n % 2 == 0, 'n must be even'

        self.n = n
        self.r = r
        self.q = q

        for i in vecpairs.keys():
            assert isinstance(i, basis_element), 'keys must be basis elements'
            assert i.n == n, 'basis elements must have n nodes'
            assert i.r == r, 'basis elements must have r anchors'
        self.vecpairs = vecpairs

    # returns coefficient of a basis element
    def coeff(self, i):
        if i in self.vecpairs.keys():
            return self.vecpairs[i]
        else:
            return 0

    # Action by generator (1+T_i), returns new element rather than editing current
    # note the use of __add__ implemented below
    def mult_generator(self, i):
        assert i > 0 and i < self.n + self.r, 'transposition out of range'
        new = module_element(self.n, {}, self.q, self.r)
        for basis in self.vecpairs.keys():
            if basis.pair(i) == i+1:
                new = new + \
                    module_element(
                        self.n, {basis: self.vecpairs[basis]*(self.q + 1)}, self.q, self.r)
            elif basis.pair(i) != i or basis.pair(i+1) != i + 1:
                match = basis.get_matches()

                if basis.pair(i) == i or basis.pair(i+1) == i+1:
                    a = match[i if basis.pair(i+1) == i+1 else i+1]
                    match[a] = a
                else:
                    a = match[i]
                    b = match[i+1]
                    match[a] = b
                    match[b] = a
                
                match[i] = i+1
                match[i+1] = i

                newmatch = {}
                for j in match.keys():
                    if match[j] != j:
                        newmatch[j] = match[j]
                newbasis = basis_element(self.n, newmatch, self.r)
                new = new + \
                    module_element(
                        self.n, {newbasis: self.vecpairs[basis]*sqrt(self.q)}, self.q,self.r)
        return new

    # Action by constant c, returns new element
    def scale(self, c):
        vecnew = {}
        for i in self.vecpairs.keys():
            vecnew[i] = self.vecpairs[i]*c
        return module_element(self.n, vecnew, self.q,self.r)

    # lets you add two elements
    def __add__(self, elt):
        vecnew = {}
        for i in self.vecpairs.keys():
            vecnew[i] = self.vecpairs[i]+elt.coeff(i)
        for i in elt.vecpairs.keys():
            if i not in self.vecpairs.keys():
                vecnew[i] = elt.coeff(i)
        return module_element(self.n, vecnew, self.q, self.r)

    # deletes all elements that are nearly zero
    # useful if you have many elements whose coefficients are zero, or
    # close because of float point error
    def round(self):
        deletions = []
        for i in self.vecpairs.keys():
            if self.vecpairs[i] < 0.0001 and self.vecpairs[i] > -0.0001:
                deletions.append(i)
        for i in deletions:
            del self.vecpairs[i]

    # to check if the element is zero
    def zero(self):
        nothing = True
        for i in self.vecpairs.keys():
            if self.vecpairs[i] > .0001 or self.vecpairs[i] < -0.0001:
                nothing = False
        return nothing

    def get_vecpairs(self):
        return self.vecpairs

    def output_latex(self, spacer="\\hspace{5pt}", delim="\\\\ \n +"):
        """
        Gives a latex command to call the Matchings macro corresponding to this element

        ARGUMENTS:
        spacer -- string placed between the coefficient and matching
        delim -- delimiter between each coefficient matching pair
        """
        matchings = ["{}{}{}".format(
            self.vecpairs[key], spacer, key.output_latex()) for key in self.vecpairs.keys()]
        return delim.join(matchings)

    def __str__(self):
        return "coeffs: " + str(self.vecpairs.values()) + " number of coeffs: "+str(len(self.vecpairs.keys()))

    def __eq__(self, i):
        return self.vecpairs == i.vecpairs

# Returns list of the C_n basis elements

def generate_basis(n, q=1, efficient_dic={}):
    if n == 2:
        return [basis_element(2, {1: 2, 2: 1})]
    if n == 0:
        return [basis_element(0, {})]
    elements = []
    for i in range(1, n//2+1):
        if (2*i-2) in efficient_dic.keys():
            sub_elts1 = efficient_dic[2*i-2]
        else:
            sub_elts1 = generate_basis(2*i-2, efficient_dic)
        if (n-2*i) in efficient_dic.keys():
            sub_elts2 = efficient_dic[n-2*i]
        else:
            sub_elts2 = generate_basis(n-2*i, efficient_dic)
        for j in sub_elts1:
            for k in sub_elts2:
                match = {1: 2*i, 2*i: 1}
                for ind in range(2*i-2):
                    match[ind+2] = j.pair(ind+1)+1
                for ind in range(n-2*i):
                    match[ind+2*i+1] = k.pair(ind+1)+2*i
                elements.append(basis_element(n, match))
    efficient_dic[n] = elements
    return elements

def generate_basis_generalized(n,q=1,r=0):
    prebasis = generate_basis(n+2*r,q)
    basis = []
    for elt in prebasis:
        quot = False
        for i in range(n+r+1,n + 2*r + 1):
            quot = quot or elt.pair(i) >= i
        if not quot:
            matches = elt.get_matches()
            newmatches = {}
            for i in range(1,n + r + 1):
                if matches[i] <= n + r:
                    newmatches[i] = matches[i]
                    newmatches[matches[i]] = i
            basis += [basis_element(n,newmatches,r)]

    return basis

