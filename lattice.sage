from sage.modules.free_module_element import vector
from time import process_time


def genXi(m,rho=12,eta=200,gam=500):
    p= 706549229#random_prime(2^eta)
    ciphertexts =  [56405845507494530020941008480572940286181689237258854,
                    39904821464460948494700284192336525523357407545067668,
                    56294991345433284900612805613249060787237279328022519]#[ p*ZZ.random_element(2^(gam-eta))+ZZ.random_element(2^rho) for _ in range(m)]

    return p,ciphertexts


# Step 1: Define the lattice L with the basis matrix
def create_lattice(ciphertexts):
    n = len(ciphertexts)
    L = matrix(n)
    for i in range(n):
        row = [int(i ==j) for j in range(n)]
        row[-1] = ciphertexts[i]
        L[i] = vector(row) 
    return L

# Step 2: Apply the LLL algorithm to reduce the basis matrix
def apply_LLL(L):
    B = L.LLL() # Apply LLL algorithm to L
    return B

def create_lattice_prime(ai,c,C=1):
    n = len(ai)
    L = matrix(n+1)
    for i in range(n+1):
        row = [int(i ==j) for j in range(n+1)]
        if i == n:
            row[-1] = C*sum([c[i] * ai[i] for i in range(n)])

            L[i] = vector(row)
            return L

        row[-1] = C*ai[i]
        L[i] = vector(row)

def attack(ciphertexts,m):
    L = create_lattice(ciphertexts)
    B = apply_LLL(L)
    ci = ciphertexts
    p=0
    for row in B:
        # Compute a_m
        am = (row[-1] - sum([row[i] * ci[i] for i in range(m-1)])) / ci[-1]
        
        row[-1] = am
        # Compute alpha
        alpha = max([log(abs(row[i]),2) for i in range(m-1)])
        # Compute rho
        rho = (1/m) * (log(abs(ci[m-1]),2) + log(sqrt(m)/(m+1),2)) - alpha
        C = int(m^((m-1)/2) * (m+1)^(-(m+1)/2) * 2^(m*rho-alpha) + 1)
        print(C)
        L = create_lattice_prime(row,ci,C=C)
        B = apply_LLL(L)
        for row in B:
            if row[-1] == 0:
                p = gcd(ci[0] - abs(row[0]),ci[1] - abs(row[1]))
                if p > 1:
                    q = [(ci[i] - row[i])/p for i in range(m)]
            return p            

def main(eta,rho,gama,m):
    p,c = genXi(m,rho,eta,gama)
    p_recover = attack(c,m)

    return p == p_recover
etas = [200]
ms = [3]
results = {}
iterations = 5
for i in range(1):
    print(i)
    for eta in etas:
        for index,m in enumerate(ms,start=1):
            cnt = 0
            t1_start = process_time() 
            gamma = eta*(2^index)
            if m == 12:
                gamma = eta*11

            for i in range(1):
                rho = 12 #int(sqrt(eta))
                b = main(eta,rho,gamma,m)
                if(b):
                    cnt += 1
            t1_stop = process_time()
            time = (t1_stop-t1_start)/100
            print(f"Success rate for eta={eta}, rho={rho}, gamma={gamma}, m={m}: {cnt}%")
            print(f"Time for 100 iterations {t1_stop-t1_start} or {time} per iteration\n")
            if (eta,m) not in results:
                results[(eta,m)] = [[cnt],[time]]
            else:
                results[(eta,m)][0].append(cnt)
                results[(eta,m)][1].append(time)
            
for k,v in results.items():
    print("params: ",k)
    print("average success rate: ", (sum(v[0])//iterations))
    print("average time: ", (sum(v[1])/iterations))
    print("\n")
    