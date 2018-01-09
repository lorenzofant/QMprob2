import matplotlib.pyplot as plt
import numpy as np
import itertools

class System():
    def __init__(self):
        self.a = 1.
        self.n = 16.
        self.epsilon = 120. # deviating from the value given in the Q 99.55 for better fitting
        self.sigma = 3.4

    def potential(self, d, vtype = "lj"):
        vval=0.
        r = d/self.sigma # to avoid recalculation of this term when calling pow()
        if vtype=="lj": vval= (pow(r, -12) - pow(r, -6))
        if vtype=="elj": vval= (pow(r, -14)+ pow(r, -12)+ pow(r, -10)+ pow(r, -8)- pow(r, -6))
        if vtype=="exp6": vval= np.exp(-1*d) - pow(r, 6)
        return self.epsilon*vval

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
            if (i, j, k) != (0, 0, 0): self.energy += self.fcccell(self.a*i, self.a*j, self.a*k, potentialType)

    def latticebcc(self, potentialType):
        self.closen = 0.
        self.energy = self.bcccell0(potentialType)
        r = np.linspace(-self.n, self.n, 2.*self.n+1)
        for i, j, k in itertools.product(r,r,r):
            if (i, j, k) != (0, 0, 0): self.energy += self.bcccell(self.a*i, self.a*j, self.a*k, potentialType)

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
        return self.potential(self.distance(a*x+(y+z)/a, y/a, z/a), potentialType)

    def bcccell(self, x, y, z, potentialType):
        a = 1./(3.**0.5)
        return self.potential(self.distance(2.*a*x, 2.*a*y, 2.*a*z)) + self.potential(self.distance(2.*a*x+self.a*a, 2.*a*y+self.a*a, 2.*a*z+self.a*a), potentialType)

    def hcpcell(self, x, y, z, potentialType):
        return self.potential(self.distance(x + y/2., y/2.*(3.**0.5), 24.**0.5*z/3.), potentialType) + self.potential(self.distance(x+(y+self.a)/2., (y+self.a/3.)/2.*(3.**0.5), 24.**0.5*(z+self.a/2.)/3.), potentialType)

    def sccell0(self, potentialType):
        return 0.

    def fcccell0(self, potentialType):
        return 0.

    def bcccell0(self, potentialType):
        a = 1./2.**0.5
        return self.potential(self.distance(self.a*a, self.a*a, self.a*a), potentialType)

    def hcpcell0(self, potentialType):
        return self.potential(self.distance(self.a/2., self.a/3./2.*(3.**0.5), 6.**0.5*self.a/3.), potentialType)


# main body of the program
a = System()
esc = []
efcc = []
ebcc = []
ehcp = []
x = np.linspace(.95, 2, 40)
for d in x:
    print( d )
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
    # print( efcc[-1] )
    a.latticebcc("lj")
    ebcc.append(a.energy)
    ehcp.append(eh/1.)
    # print( ehcp[-1] )
# plt.plot(x, esc, label='SC')
# plt.plot(x, efcc, label='FCC')
# print( sum(efcc)/len(efcc) )
# plt.plot(x, ebcc, label='BCC')
# print( sum(ehcp)/len(ehcp) )
# plt.plot(x, ehcp, label='HCP')
# plt.grid()
# plt.legend()
# plt.xlabel('Lattice spacing [$\sigma$]')
# plt.ylabel('Energy [$\epsilon$]')
# plt.show()
# np.savetxt('energie6', [x, efcc, ehcp])
print("sc min = ", min(esc), "bcc min = ", min(ebcc), ":: fcc min = ", min(efcc), ":: hcp min = ", min(ehcp))
