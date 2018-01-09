import matplotlib.pyplot as plt
import numpy as np
import itertools

class System():
    def __init__(self):
        self.a = 1.
        self.n = 16.
        self.epsilon = 1.
        self.sigma = 1.

    def potential(self, d, vtype = "lj"):
        vval=0.
        if vtype=="lj": vval= 4.*self.epsilon*(pow(self.sigma/d, 12) - pow(self.sigma/d, 6))
        if vtype=="elj": vval= 4.*self.epsilon*(pow(self.sigma/d, 14)
                                                - pow(self.sigma/d, 13)
                                                + pow(self.sigma/d, 12)
                                                - pow(self.sigma/d, 11)
                                                + pow(self.sigma/d, 10)
                                                - pow(self.sigma/d, 9)
                                                + pow(self.sigma/d, 8)
                                                - pow(self.sigma/d, 6))
        if vtype=="exp6": vval= self.A1*np.exp(-1*self.B1*d) - pow(self.sigma/d, 6)
        return vval

    def distance(self, x, y, z):
        d = (x**2 + y**2 + z**2)**0.5
        if d < self.a*1.1:
            self.closen += 1
        return d

    def latticesc(self, potentialType):
        self.closen = 0.
        self.energy = self.sccell0(potentialType)
        r = np.linspace(-self.n, self.n, 2.*self.n+1)
        for i, j, k in itertools.product(r,r,r):
            if (i,j,k) != (0,0,0): self.energy += self.sccell(self.a*i, self.a*j, self.a*k, potentialType)

    def latticefcc(self, potentialType):
        self.closen = 0.
        #self.a = 1.5416
        self.energy = self.fcccell0(potentialType)
        r = np.linspace(-self.n, self.n, 2.*self.n+1)
        for i, j, k in itertools.product(r,r,r):
            if (i,j,k) != (0,0,0): self.energy += self.fcccell(self.a*i, self.a*j, self.a*k, potentialType)

    def latticebcc(self, potentialType):
        self.closen = 0.
        self.energy = self.bcccell0(potentialType)
        r = np.linspace(-self.n, self.n, 2.*self.n+1)
        for i, j, k in itertools.product(r,r,r):
            if (i,j,k) != (0,0,0): self.energy += self.bcccell(self.a*i, self.a*j, self.a*k, potentialType)

    def latticehcp(self, potentialType):
        self.closen = 0.
        #self.a = 1.09
        self.energy = self.hcpcell0(potentialType)
        r = np.linspace(-self.n, self.n, 2.*self.n+1)
        r2 = np.linspace(-self.n/2., self.n/2., self.n+1)
        for i, j, k in itertools.product(r,r,r2):
            if (i,j,k) != (0,0,0): self.energy += self.hcpcell(self.a*i, self.a*j, self.a*k, potentialType)

    def sccell(self, x, y, z, potentialType):
        return self.potential(self.distance(x, y, z), potentialType)

    def fcccell(self, x, y, z, potentialType):
        a = 2.**0.5
        return self.potential(self.distance(a*x+(y+z)/a + self.noise*np.random.normal(), y/a + self.noise*np.random.normal(), z/a + self.noise*np.random.normal()), potentialType)

    def bcccell(self, x, y, z, potentialType):
        a = 1./(3.**0.5)
        return self.potential(self.distance(2.*a*x, 2.*a*y, 2.*a*z)) + self.potential(self.distance(2.*a*x+self.a*a, 2.*a*y+self.a*a, 2.*a*z+self.a*a), potentialType)

    def hcpcell(self, x, y, z, potentialType):
        return self.potential(self.distance(x + y/2. + self.noise*np.random.normal(), y/2.*(3.**0.5)+self.noise*np.random.normal(), 24.**0.5*z/3.+self.noise*np.random.normal()), potentialType) + self.potential(self.distance(x+(y+self.a)/2.+self.noise*np.random.normal(), (y+self.a/3.)/2.*(3.**0.5)+self.noise*np.random.normal(), 24.**0.5*(z+self.a/2.)/3.+self.noise*np.random.normal()), potentialType)

    def sccell0(self, potentialType):
        return 0.

    def fcccell0(self, potentialType):
        return 0.

    def bcccell0(self, potentialType):
        a = 1./2.**0.5
        return self.potential(self.distance(self.a*a, self.a*a, self.a*a), potentialType)

    def hcpcell0(self, potentialType):
        return self.potential(self.distance(self.a/2.+self.noise*np.random.normal(), self.a/3./2.*(3.**0.5)+self.noise*np.random.normal(), 6.**0.5*self.a/3.+self.noise*np.random.normal()), potentialType)


# main body of the program
a = System()
esc = []
efcc = []
ebcc = []
ehcp = []
x = np.linspace(.95, 1.2, 41)
for d in x:
    print( d )
    a.noise = 0.
    a.a = d
    ef = 0.
    eh = 0.
    a.latticesc("lj")
    esc.append(a.energy)
    for i in range(1):
        #print( i )
        #a.a = 1.54179
        a.latticefcc("lj")
        ef += (a.energy)
        #a.a = 1.0902
        a.latticehcp("lj")
        eh += (a.energy)
    efcc.append(ef/1.)
    print( efcc[-1] )
    a.latticebcc("lj")
    ebcc.append(a.energy)
    ehcp.append(eh/1.)
    print( ehcp[-1] )
plt.plot(x, esc, label='SC')
plt.plot(x, efcc, label='FCC')
print( sum(efcc)/len(efcc) )
plt.plot(x, ebcc, label='BCC')
print( sum(ehcp)/len(ehcp) )
plt.plot(x, ehcp, label='HCP')
plt.grid()
plt.legend()
plt.xlabel('Lattice spacing [$\sigma$]')
plt.ylabel('Energy [$\epsilon$]')
plt.show()
# np.savetxt('energie6', [x, efcc, ehcp])
# print("sc min = ", min(esc), "bcc min = ", min(ebcc), ":: fcc min = ", min(efcc), ":: hcp min = ", min(ehcp))
