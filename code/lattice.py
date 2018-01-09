import matplotlib.pyplot as plt
import numpy as np
import itertools

class System():
    def __init__(self):
        self.a = 3.75
        self.n = 16.
        self.epsilon = 99.55 # deviating from the value given in the Q 99.55 for better fitting
        self.sigma = 3.4
        # constants for BFW potential
        self.bfw_epsilon=142.1; self.bfw_sigma=3.76
        self.bfw_C6 = -1.10727; self.bfw_A0 = 0.27783; self.bfw_A3 = -25.2696
        self.bfw_C8 = -0.16971; self.bfw_A1 = -4.50431; self.bfw_A4 = -102.0195
        self.bfw_alpha = 12.5; self.bfw_C10 = -0.01361; self.bfw_A2 = -8.33122; self.bfw_A5 = -113.25
        self.bfw_delta = 0.01

    def potential(self, d, vtype = "lj"):
        vval=0.
        if vtype=="lj": 
            r = d/self.sigma # to avoid recalculation of this term when calling pow()
            vval= 4*self.epsilon*(pow(r, -12) - pow(r, -6))
        if vtype=="bfw": 
            r = d/self.bfw_sigma
            vval= self.bfw_epsilon*(np.exp(self.bfw_alpha*(1-r))*(self.bfw_A1*(r-1)+self.bfw_A2*pow((r-1),2)+self.bfw_A3*pow((r-1),3)+self.bfw_A4*pow((r-1),4)+self.bfw_A5*pow((r-1),5))+(self.bfw_C6/(self.bfw_delta+pow(r,6))+self.bfw_C8/(self.bfw_delta+pow(r,8))+self.bfw_C10/(self.bfw_delta+pow(r,10))))
        return vval

    def distance(self, x, y, z):
        d = (x**2 + y**2 + z**2)**0.5
        if d < self.a*1.1:
            self.closen += 1
        return d

    def lattice(self, latticeType, potentialType):
        self.closen = 0.
        r = np.linspace(-self.n, self.n, 2.*self.n+1)
        r2 = np.linspace(-self.n/2., self.n/2., self.n+1)
        if latticeType=="sc":
            self.energy = self.sccell0(potentialType)
            for i, j, k in itertools.product(r,r,r):
                if (i,j,k) != (0,0,0): self.energy += self.sccell(self.a*i, self.a*j, self.a*k, potentialType)
        if latticeType == "fcc":
            self.energy = self.fcccell0(potentialType)
            for i, j, k in itertools.product(r,r,r):
                if (i, j, k) != (0, 0, 0): self.energy += self.fcccell(self.a*i, self.a*j, self.a*k, potentialType)
        if latticeType == "bcc":
            self.energy = self.bcccell0(potentialType)
            for i, j, k in itertools.product(r,r,r):
                if (i, j, k) != (0, 0, 0): self.energy += self.bcccell(self.a*i, self.a*j, self.a*k, potentialType)
        if latticeType == "hcp":
            self.energy = self.hcpcell0(potentialType)
            for i, j, k in itertools.product(r,r,r2):
                if (i,j,k) != (0,0,0): self.energy += self.hcpcell(self.a*i, self.a*j, self.a*k, potentialType)

    def sccell(self, x, y, z, potentialType):
        a=1.0
        return self.potential(self.distance(a*x, a*y, a*z), potentialType)

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
x = np.linspace(3, 5, 40)
for d in x:
    print( d )
    a.a = d
    ef = 0.
    eh = 0.
    a.lattice("sc", "lj")
    esc.append(a.energy)
    for i in range(1):
        #print( i )
        #a.a = 1.54179
        a.lattice("fcc", "lj")
        ef += (a.energy)
        #a.a = 1.0902
        a.lattice("hcp", "lj")
        eh += (a.energy)
    efcc.append(ef/1.)
    # print( efcc[-1] )
    a.lattice("bcc", "lj")
    ebcc.append(a.energy)
    ehcp.append(eh/1.)
    # print( ehcp[-1] )
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
print("sc min = ", min(esc), "bcc min = ", min(ebcc), ":: fcc min = ", min(efcc), ":: hcp min = ", min(ehcp))
