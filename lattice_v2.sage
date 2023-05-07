from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()from sage.modules.free_module_element import vector

lambd = 2^100
def genXi(rho=10,eta=100,gam=200):
    x1 = 5
    x0 = 1
    while x1 > x0 or gcd(x0,x1) !=1:
        p=random_prime(2^eta)
        r= ZZ.random_element(2^rho) 
        q0=ZZ.random_element(2^(gam-eta-1),2^(gam-eta))
        q1=ZZ.random_element(2^(gam-eta))

        x0 = p*q0
        x1 = p * q1 + r


    return x0,x1,p

def enc(pk,m,n=32):
    d=2*n
    N=ZZ.random_element(2^(d-2),2^d) * 2
 
    c = (m + N * pk[1]) % pk[0]
    return c

# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    # get the size of the matrix
    t = len(ciphertexts)
    # create the (t, 2) matrix with a constant and the vector u
    A = matrix([[lambd*ciphertexts[i]] for i in range(t)])

    # create the identity matrix of size (t, t-1)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M

def create_lattice_prime(t,u):
    # create the matrix of size (t, t-3)
    A = matrix([[lambd * u[i][j] for i in range(t-3)] for j in range(t)])

    # create the identity matrix of size (t, t)
    I = identity_matrix(t)

    # concatenate the two matrices horizontally
    M = A.augment(I)
    return M


def attack():
    t = 5
    x0,x1,p=genXi()
    m = [ZZ.random_element(2) for i in range(t)]
    ciphertexts = [enc((x0,x1),mi) for mi in m]
    print(m)
    L = create_lattice(ciphertexts)
    B = L.LLL()
    B = B.delete_columns([0])
    B = B.delete_rows([B.nrows()-1])
    M = create_lattice_prime(t,B)
    B = M.LLL()
    delete_columns = B.ncols() - B.nrows()
    for i in range(delete_columns):
        B = B.delete_columns([0])

    print(B)

def main():
   attack()

main()