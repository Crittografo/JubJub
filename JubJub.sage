# twisted Edwards curve has form  
# e*x^2 + y^2 = 1 + d*x^2*y^2 mod p


# BabyJubJub is Edwards curve related to curve Alt_Bn128
# This curve is optimal to be used in zkSnarks in Ethereum

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


# JubJub is Edwards curve related to BLS12-381 curves (that used in Zcash Sapling) and has modulus p equal to prime subgroup order of BLS12-381 
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


def add_points(x1, y1, x2, y2, d, e, F):
    xx = F(x1*x2)
    yy1 = F(y1)*F(y2)
    dxxyy = F(d*xx*yy1)
    
    x3 = (F(x1*y2) + F(y1*x2)) * (1 + dxxyy)^-1
    
    inv2 = (F(1) - dxxyy)^-1 
    exx = F(e)*xx
        
    y3 = (yy1 - exx) * inv2
    
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

def test_edw_addition(d, e, p):
    
    Fp = GF(p)
    
    a, b = coeff_edw_to_w(d, e, Fp)
    
    E = EllipticCurve(Fp, [a, b]) 
    
    R1 = E.random_point()
    R2 = E.random_point()
    R3 = R1 + R2
    
    r1x, r1y = coord_w_to_edw(R1[0], R1[1], d, e, Fp)
    r2x, r2y = coord_w_to_edw(R2[0], R2[1], d, e, Fp)
    r3x, r3y = coord_w_to_edw(R3[0], R3[1], d, e, Fp)
        
    r3x_4check, r3y_4check = add_points(r1x, r1y, r2x, r2y, d, e, Fp)
    
    if r3x == r3x_4check and r3y == r3y_4check :
        return True
    else:
        return False
    
        
    
def test_sign_verify(testN, p, q, d, e, Gx, Gy):

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
                            
            points.append([F(x), F(y)])
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


def get_precomp_point(points, window_value, bits_in_window, d, e, F):
    x = 0
    y = 1

    for i in range(bits_in_window):
        if (1 << i) & window_value:
            x, y = add_points(x, y, points[i][0], points[i][1], d, e, F)
            
    return x, y

# input is one of chunks chunk_points: all generators for one of chunks
# returns array with all possible results of addition of generators, where inde[ in array is value of window 

def get_points_for_all_window_values(chunk_generators, bits_in_window, d, e, F):
    
    # values of chunk for all possible values of window (from 0 to (2^bits_in_window)-1 )
    chunk_values_of_points = []
      
    for window_value in range(2^bits_in_window):
   
        x, y = get_precomp_point(chunk_generators, window_value, bits_in_window, d, e, F)
        chunk_values_of_points.append([x, y])
              
    return chunk_values_of_points

def extractKBits(num,k,p): 
  
     # convert number into binary first 
     binary = bin(num) 
  
     # remove first two characters 
     binary = binary[2:] 
     
     zeros2add = k - len(binary)%k
     
     while zeros2add != 0 :
          binary = "0" + binary
          zeros2add -=1

     end = len(binary) - p 
     start = end - k + 1
  
     # extract k  bit sub-string 
     kBitSubStr = binary[start : end+1] 
  
     # convert extracted sub-string into decimal again 
     #(int(kBitSubStr,2)) 

     return int(kBitSubStr,2)


def pedersen_hash_ex(data, window_len, chanks_precomps, d, e, F):
    x = 0
    y = 1
    
    number = 0    
    b = bytearray(data)
    pow = 2^8
    lenth = len(b)
    
    for i in range(lenth):
        number += b[i]*(pow^i)
        
    length_in_bits = number.nbits()
    
    number_of_chunks = (length_in_bits)//window_len
    
    if (length_in_bits%window_len) != 0:
        number_of_chunks += 1 
     
    for j in range(number_of_chunks):
        Kbit = extractKBits(number, window_len, 1 + j*window_len)
        precompX = chanks_precomps[j][Kbit][0]
        precompY = chanks_precomps[j][Kbit][1]
        x, y = add_points(x, y, precompX, precompY, d, e, F)
        #print hex(Kbit)
    
    return x, y


#print test_sign_verify(5, p_BabyJubJub, q_BabyJubJub, d_BabyJubJub, e_BabyJubJub,Gx_BabyJubJub, Gy_BabyJubJub )
#print test_sign_verify(5, p_JubJub, q_JubJub, d_JubJub, e_JubJub,Gx_JubJub, Gy_JubJub )
#test_edw_addition(d_BabyJubJub, e_BabyJubJub, p_BabyJubJub)


# test pedersen hash: 

seed = "Hello, World!"

# number of generators must be divisible by window len
window_len = 8  

generators_list = create_points_for_hash(seed, 512 + window_len, d_BabyJubJub, e_BabyJubJub, F_BabyJubJub, p_BabyJubJub)
precomp_table_of_chunks = []

num_of_chunks = len(generators_list)//window_len

for i in range(num_of_chunks):
    
    ith_chunk_precomp = get_points_for_all_window_values(generators_list[i*window_len:], window_len, d_BabyJubJub, e_BabyJubJub, F_BabyJubJub)
    precomp_table_of_chunks.append(ith_chunk_precomp)
    
print "OK. chunks calculated"

data = [0x34, 0x11, 0x1c]

x, y = pedersen_hash_ex(data, window_len, precomp_table_of_chunks, d_BabyJubJub, e_BabyJubJub, F_BabyJubJub)
print(x, y)

x1, y1 = pedersen_hash(data, generators_list, d_BabyJubJub, e_BabyJubJub, F_BabyJubJub)
print(x1, y1) 

if x == x1 and y == y1:
    print "test OK"
else:
    print "test failed"
