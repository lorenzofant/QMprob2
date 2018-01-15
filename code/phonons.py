import matplotlib.pyplot as plt
import numpy as np
import itertools

class System():
    def __init__(self):
        self.n = 2.
        # constants for LJ potential
        self.epsilon = 99.55 #137.37 #
        self.sigma = 3.758/(2.**(1./6.))
        self.a = 1.0902*self.sigma #3.75
        self.mass = 1.
        
    
    def potentialp(self, x, i, j):
        d = np.linalg.norm(x)
        return self.LJpotential_deriv2(d)*(np.vdot(x,i)*np.vdot(x,j))/(d**2) #+ self.LJpotential_deriv(d)*(np.vdot(i,j)/d - (np.vdot(x,i)*np.vdot(x,j))/(d**3))

    def LJpotential_deriv(self, d):
        r = d/self.sigma
        return -4.*self.epsilon/self.sigma*(12.*pow(r, -13) - 6.*pow(r, -7))

    def LJpotential_deriv2(self, d):
        r = d/self.sigma
        return 4.*self.epsilon/(self.sigma**2)*(156.*pow(r, -14) - 42.*pow(r, -8))

    def dynamicmatrix(self):
        self.dm = []
        r = np.linspace(-self.n, self.n, 2.*self.n+1)
        r2 = np.linspace(-self.n/2., self.n/2., self.n+1)
        for l in range(3):
            c = np.array([0.,0.,0.])
            c[l]+=1.
            a = []
            for m in range(3):
                d = np.array([0.,0.,0.])
                d[m]+=1.
                self.energy = 0.#self.hcpcell0(potentialType)
                for i, j, k in itertools.product(r,r,r):
                    x = self.a*np.array([i, j, k])
                    if (i,j,k) != (0.,0.,0.): 
                        self.energy += self.fcccellp(x, c, d)*(1.-self.phasek(self.fccvector(x)))
                self.energy /= self.mass
                a.append(self.energy)
            self.dm.append(a)
        self.w = np.sort(np.linalg.eigvals(self.dm))

    def dynamicmatrixhcp(self):
        self.dm = []
        r = np.linspace(-self.n, self.n, 2.*self.n+1.)
        r2 = np.linspace(-self.n/2., self.n/2., self.n+1.)
        for l in range(6):
            c = np.array([0.,0.,0.])
            c[l%3]+=1.
            a = []
            for m in range(6):
                d = np.array([0.,0.,0.])
                d[m%3]+=1.
                self.energy = 0.#self.hcpcell0(potentialType)
                for i, j, k in itertools.product(r,r,r2):
                    x = self.a*np.array([i, j, k])
                    if l//3.==m//3.:
                        if (i,j,k) != (0,0,0):              
                             self.energy += self.hcpcellp1(x, c, d, -l//3.+m//3.)*(2.*(1.-l//3.)*(1.-m//3.)+2.*l//3.*m//3.-self.phasek(self.hcpvector(x,-l//3.+m//3.)))
                    else:
                        self.energy += self.hcpcellp1(x, c, d, -l//3.+m//3.)*(2.*(1.-l//3.)*(1.-m//3.)+2.*l//3.*m//3.-self.phasek(self.hcpvector(x,-l//3.+m//3.)))
                self.energy /= self.mass
                a.append(self.energy)
            self.dm.append(a)
        print self.dm
        print(np.shape(self.dm))
        self.w = np.sort(np.linalg.eigvals(np.transpose(self.dm)))

    def fccvector(self,x):
        a = 2.**0.5
        return np.array([(x[0]+x[1])/a, (x[0]+x[2])/a, (x[1]+x[2])/a])

    def hcpvector(self,x,e):
        return[x[0] + (x[1] + e*self.a)/2., (x[1] + e*self.a/3.)/2.*(3.**0.5), 24.**0.5*(x[2] + e*self.a/2.)/3.]
    def phasek(self, x):
        return np.e**(1.j*np.dot(self.k,x))#np.cos(np.dot(self.k,x))

    def distance(self, x, y, z):
        d = np.linalg.norm(x)
        return d


    def hcpcellp2(self, x, c, d):
        return self.potentialp([x[0] + x[1]/2., x[1]/2.*(3.**0.5), 24.**0.5*x[2]/3.], c, d) 

    def hcpcellp1(self, x, c, d, e):
        return self.potentialp(self.hcpvector(x,e), c, d)

    def fcccellp(self, x, c, d):
        return self.potentialp(self.fccvector(x), c, d)


# main body of the program
a = System()
esc = []; ebcc = []
efcc_lj = []; efcc_bfw = []; efcc_hfd = []; efcc_bbms = []
ehcp_lj = []; ehcp_bfw = []; ehcp_hfd = []; ehcp_bbms = []
x = np.linspace(0.,2.**0.5/a.a*np.pi,100) 
dx = x[1]-x[0]
w = []
for el in x:
    print el
    a.k = el*np.array([0.,1.,0.])#np.array([0.5,-0.5,-0.5])
    a.dynamicmatrix()
    w.append(a.w)
for el in x[:50]:
    print el
    a.k = x[-1]*np.array([el/x[-1],1.,0.])#np.array([0.5,-0.5,-0.5])
    a.dynamicmatrix()
    w.append(a.w)
for el in x[:25]:
    print el
    a.k = x[-1]*np.array([0.5+el/x[-1],1.-el/x[-1],0.])#np.array([0.5,-0.5,-0.5])
    a.dynamicmatrix()
    w.append(a.w)
for el in x[:100][::-1]:
    print el
    a.k = el*np.array([.75,.75,0.])#np.array([0.5,-0.5,-0.5])
    a.dynamicmatrix()
    w.append(a.w)
w = np.array(np.transpose(w))
y = np.array(range(275))*2.**0.5/a.a*np.pi/100.
print w
for el in w:
    plt.plot(y,el**0.5)
    print (el[1]-el[0])/dx
    print (el[-1]-el[-2])/dx
plt.grid()
plt.xlabel('k')
plt.ylabel('$\omega$')
plt.show()
