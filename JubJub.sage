# twisted Edwards curve has form  
# e*x^2 + y^2 = 1 + d*x^2*y^2 mod p

# JubJub is Edwards curve. 
# JubJub curve is related to BLS12-381 curves and has modulus p equal to prime subgroup order of BLS12-381  
# that helps to create optimal zkSnarks with high (128-bit) level of security and minimal number of constraints 
# for any cryptosystem whose operations are in field GF(p), where p is a modulus of JubJub   
# e = -1, d = -10240/10241 mod p  ( == 0x2a9318e74bfa2b48f5fd9207e6bd7fd4292d7f6d37579d2601065fd6d6343eb1)
# p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001 
# subgroup order: 
q_JubJub = 6554484396890773809930967563523245729705921265872317281365359162392183254199
d_JubJub = 0x2a9318e74bfa2b48f5fd9207e6bd7fd4292d7f6d37579d2601065fd6d6343eb1
e_JubJub = -1
p_JubJub = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
#p_JubJub = 52435875175126190479447740508185965837690552500527637822603658699938581184513
    
F_JubJub = GF(p_JubJub)
Gx_JubJub = F_JubJub(8076246640662884909881801758704306714034609987455869804520522091855516602923)
Gy_JubJub = F_JubJub(13262374693698910701929044844600465831413122818447359594527400194675274060458)


# BabyJubJub is Edwards curve with parameters 
# BabyJubJub curve is related to BN128
# subgroup order:
q_BabyJubJub = 2736030358979909402780800718157159386076813972158567259200215660948447373041
	      
#coefficients e and d: 
d_BabyJubJub = 168696 
e_BabyJubJub = 168700
#modulus
p_BabyJubJub = 21888242871839275222246405745257275088548364400416034343698204186575808495617 
F_BabyJubJub = GF(p_BabyJubJub)

Gx_BabyJubJub = F_BabyJubJub(17777552123799933955779906779655732241715742912184938656739573121738514868268)
Gy_BabyJubJub = F_BabyJubJub(2626589144620713026669568689430873010625803728049924121243784502389097019475)
 
 
def add_points(x1, y1, x2, y2, d, e, F):
    xx = F(x1*x2)
    yy = F(y1*y2)
    dxxyy = F(d*xx*yy)
    
    x3 = (F(x1*y2) + F(y1*x2)) * (1 + dxxyy)^-1
    y3 = (xx - F(e)*yy) * (1 - dxxyy)^-1
    return x3, y3; 


def calc_st(d, e, F):
    
    s = (F(e) - F(d))*F(4)^-1 
    t = (F(e) + F(d))*F(6)^-1 
    
    return s, t
 
  

def coeff_edw_to_w(d, e, F):
    
    s, t = calc_st(d, e, F)
    
    a = s^2 - 3*t^2
    b = 2*t^3 - t*s^2
    
    return a, b


    
def coord_edw_to_w(u, v, d, e, F):   
     
    s, t = calc_st(d, e, F)
    
    x = s*(1 + F(v))*F(1 - v)^-1 + t
    y = s*(1 + F(v))*F((1 - v)*u)^-1
       
    return x, y



def coord_w_to_edw(x, y, d, e, F):  
      
    s, t = calc_st(d, e, F)
    
    u = (F(x) - t)*F(y)^-1
    v = (F(x) - t - s)*F(x - t + s)^-1
    
    return u, v 



def paramset_edw_to_w(q, Gx, Gy, d, e, F):

    Gx_w, Gy_w = coord_edw_to_w(Gx, Gy, d, e, F)
    a, b = coeff_edw_to_w(d, e, F)

    return Gx_w, Gy_w, a, b



def scalarMultEdw(x, Px, Py, q, d, e, F ):

     # Q = x*P : 
     Fq = GF(q)
     Px_w, Py_w, a, b = paramset_edw_to_w(q, Px, Py, d, e, F)
     E = EllipticCurve(F, [a, b]) 
     Q = int(x)*E(Px_w, Py_w)
     QXEdw, QYEdw = coord_w_to_edw(Q[0], Q[1], d, e, F)

     return QXEdw, QYEdw



def gen_keypairEdv(q, Gx, Gy, d, e, F):

    Fq = GF(q)
    
    priv = Fq.random_element()
    
    PubXEdw, PubYEdw = scalarMultEdw(priv, Gx, Gy, q, d, e, F)
    
    return priv, PubXEdw, PubYEdw



def signEdw(hash, priv, q, Gx, Gy, d, e, p):

    Fq = GF(q)
    Fp = GF(p)
    recID = 0     
    
    k, rx, ry = gen_keypairEdv(q, Gx, Gy, d, e, Fp)
    s = (k^-1)*(Fq(hash) + priv*Fq(rx))
    
    recID = (int(rx))//q

   
    if ry > p//2:
         recID = recID + 0x80
   
 
    return Fq(rx), Fq(s), recID 


def verifyEdw(hash, r, s, Pubx, Puby, q, Gx, Gy, d, e, p):
    Fq = GF(q)
    Fp = GF(p)
    
    if Fq(hash) == 0:
        return False
 
    if  Fq(r) == 0:
        return False
        
    if  Fq(s) == 0:
        return False    
         
    u1 = Fq(hash)*(Fq(s)^-1)
    u2 = Fq(r)*(Fq(s)^-1)
    
    U1x, U1y = scalarMultEdw(u1, Gx, Gy, q, d, e, Fp)
    U2x, U2y = scalarMultEdw(u2, Pubx, Puby, q, d, e, Fp)
    
    Ux, Uy = add_points(U1x, U1y, U2x, U2y, d, e, Fp)
    
    if(Fq(Ux) == r):
        return True
    else:
        return False

def recover_pub_key(hash, r, s, recoveryID, Gx, Gy, q, d, e, p):
    Fq = GF(q)
    Fp = GF(p)

    # recover x:
    Rx = r + (recoveryID & 0b1111111)*q

    # recover y:
    squareRy = (1 - Fp(e)*Fp(Rx)^2)*( 1 - Fp(d)*Fp(Rx)^2)^-1
        
    if kronecker(squareRy, p) == 1:
        Ry = squareRy.sqrt()
        
        if recoveryID & 0x80:
            Ry = p - Ry
    else:
        return ValueError("Square root doesn't exist") 
            
    #Q = r^-1*s*R - r^-1*hash*G        
    
    s1 = int( Fq(r)^-1*Fq(s) )
    
    #scalarMultEdw(x, Px, Py, q, d, e, F )
    
    Q1x, Q1y = scalarMultEdw(s1, Rx, Ry, q, d, e, Fp)
    
    s2 = int( Fq(r)^-1*Fq(hash) )
    s2 = q - s2
    
    Q2x, Q2y = scalarMultEdw(s2, Gx, Gy, q, d, e, Fp)
    
    x, y = add_points(Q1x, Q1y, Q2x, Q2y, d, e, Fp)
            
    return x, y
    
    
def test(testN, p, q, e, d, Gx, Gy):

    Fp = GF(p)
    Fq = GF(q)
 
    recID = 0
    
    while testN != 0:
  
        hash = Fq.random_element()
       
        priv, pubX, pubY = gen_keypairEdv(q, Gx, Gy, d, e, Fp)
        
        u = pubX
        v = pubY
        b1 = (e*u^2 + v^2 == 1 + d*u^2*v^2)
        
        if b1 == False:
            print "gen_keypairEdv failed"
            return False 
        
        r, s, recID = signEdw(hash, priv, q, Gx, Gy, d, e, p)
  
           
        b2 = verifyEdw(hash, r, s, pubX, pubY, q, Gx, Gy, d, e, p)
        
        if b2 == False:
            print "verifyEdw failed"
            return False
        
        testN = testN - 1
        
    return True


import hashlib

def create_points_for_hash(seed, N, d, e, F, p):
    
    n = 0
     
    y = int( hashlib.sha256(seed).hexdigest(), 16)
    points = []
    
    while n < N:
        
        squareX = (F(y)^2 - 1)*(F(d)*F(y)^2 - F(e))^-1
        
        if kronecker(squareX, p) == 1:
            
            x = squareX.sqrt()
                            
            points.append([x, y])
            n = n + 1
         
        y = int( hashlib.sha256(hex(y)).hexdigest(), 16)
        
            
    return points

def pedersen_hash(data, points, d, e, F):
    x = 0
    y = 1
    
    if(len(points) == 0):
        ValueError("Array of points is empty") 
    
    if(len(points) < 8*len(data) or len(points) == 0):
        ValueError("Number of points is not enougth to calculate hash")
        
    b = bytearray(data)

    for i in range(len(b)):
        for j in range(8):
            if (b[i] >> j)&1:
                x, y = add_points(x, y, points[i*8 + j][0], points[i*8 + j][1], d, e, F)
    
    return x, y


def test_pedersen():
    seed = "Hello, World!"
    generators_list = create_points_for_hash(seed, 512, d_JubJub, e_JubJub, F_JubJub, p_JubJub)
    
    message = "The quick brown fox jumps over the lazy dog"
    x, y = pedersen_hash(message, generators_list, d_JubJub, e_JubJub, F_JubJub)
    print(x, y)   


print test(5, p_BabyJubJub, q_BabyJubJub, e_BabyJubJub, d_BabyJubJub,Gx_BabyJubJub, Gy_BabyJubJub )
print test(5, p_JubJub, q_JubJub, e_JubJub, d_JubJub,Gx_JubJub, Gy_JubJub )




