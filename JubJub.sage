# edwards curve: 
# e*x^2 + y^2 = 1 + d*x^2*y^2 mod p

# jubjub is edwards curve with parameters 
# e = -1, d = -10240/10241 mod p  ( == 0x2a9318e74bfa2b48f5fd9207e6bd7fd4292d7f6d37579d2601065fd6d6343eb1)
# p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001 
# subgroup order: 
q_JubJub = 6554484396890773809930967563523245729705921265872317281365359162392183254199
d_JubJub = 0x2a9318e74bfa2b48f5fd9207e6bd7fd4292d7f6d37579d2601065fd6d6343eb1
e_JubJub = -1
p_JubJub = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
    
F_JubJub = GF(p_JubJub)
Gx_JubJub = F_JubJub(8076246640662884909881801758704306714034609987455869804520522091855516602923)
Gy_JubJub = F_JubJub(13262374693698910701929044844600465831413122818447359594527400194675274060458)

 
 
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



def signEdw(hash, priv, q, Gx, Gy, d, e, F):
    Fq = GF(q)
         
    k, rx, ry = gen_keypairEdv(q, Gx, Gy, d, e, F)
    s = (k^-1)*(Fq(hash) + priv*Fq(rx))
    
    return Fq(rx), Fq(s)


def verifyEdw(hash, r, s, Pubx, Puby, q, Gx, Gy, d, e, F):
    Fq = GF(q)
    
    if Fq(hash) == 0:
        return False
 
    if  Fq(r) == 0:
        return False
        
    if  Fq(s) == 0:
        return False    
         
    u1 = Fq(hash)*(Fq(s)^-1)
    u2 = Fq(r)*(Fq(s)^-1)
    
    U1x, U1y = scalarMultEdw(u1, Gx, Gy, q, d, e, F)
    U2x, U2y = scalarMultEdw(u2, Pubx, Puby, q, d, e, F)
    
    Ux, Uy = add_points(U1x, U1y, U2x, U2y, d, e, F)
    
    if(Fq(Ux) == r):
        return True
    else:
        return False


def test():
    Fq = GF(q_JubJub)
    testN = 10
    
    while testN != 0:
  
        hash = Fq.random_element()
       
        priv, pubX, pubY = gen_keypairEdv(q_JubJub, Gx_JubJub, Gy_JubJub, d_JubJub, e_JubJub, F_JubJub)
        #print pubX, pubY
        u = pubX
        v = pubY
        b1 = (-u^2 + v^2 == 1 + d_JubJub*u^2*v^2)
        
        if b1 == False:
            print "gen_keypairEdv failes"
            return False 
        
        r, s = signEdw(hash, priv, q_JubJub, Gx_JubJub, Gy_JubJub, d_JubJub, e_JubJub, F_JubJub)
        
        b2 = verifyEdw(hash, r, s, pubX, pubY, q_JubJub, Gx_JubJub, Gy_JubJub, d_JubJub, e_JubJub, F_JubJub)
        
        if b2 == False:
            print "verifyEdw failed"
            return False 
        
        testN = testN - 1
        
    return True


print test()





