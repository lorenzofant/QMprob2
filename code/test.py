import numpy as np
import matplotlib.pyplot as plt
class System:
    def __init__(self):
        self.a=1.
        self.n=16.
        self.A1=1.
        self.B1=1.
        self.epsilon=1.0
        self.sigma=1.0
    
    def distance(self, x, y, z):
        return (x**2 + y**2 + z**2)**0.5
    
    def potential(self, d, vtype="lj"):
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
  
    
    def latticesc(self):
        self.energy = 0.
        for i in np.linspace(-self.n, self.n, 2.*self.n+1):
            for j in np.linspace(-self.n, self.n, 2.*self.n+1):
                for k in np.linspace(-self.n, self.n, 2.*self.n+1):
                    if (i, j, k) != (0., 0., 0.): self.energy += self.sccell(self.a*i, self.a*j, self.a*k)

    def latticefcc(self):
        self.energy = 0.
        for i in np.linspace(-self.n, self.n, 2.*self.n+1):
            for j in np.linspace(-self.n, self.n, 2.*self.n+1):
                for k in np.linspace(-self.n, self.n, 2.*self.n+1):
                    if (i,j,k) != (0,0,0): self.energy += self.fcccell(self.a*i, self.a*j, self.a*k)
                        
    def latticehcp(self):
        self.energy = self.hcpcell0()
        for i in np.linspace(-self.n, self.n, 2.*self.n+1.):
            for j in np.linspace(-self.n, self.n, 2.*self.n+1.):
                for k in np.linspace(-self.n/2., self.n/2., self.n+1.):
                    if (i, j, k) != (0, 0, 0): self.energy += self.hcpcell(self.a*i, self.a*j, self.a*k)       
                        
    def fcccell(self, x, y, z):
        a = 2.**0.5
        return self.potential(self.distance(a*x+(y+z)/a , y/a , z/a ), "lj")

    def sccell(self, x, y, z):
        return self.potential(self.distance(x, y, z),"lj")

    def hcpcell(self, x, y, z):
        vval = self.potential(self.distance(x + y/2., y/2.*(3.**0.5), 24.**0.5*z/3.), "lj") 
        + self.potential(self.distance(x+(y+self.a)/2., (y+self.a/3.)/2.*(3.**0.5), 24.**0.5*(z+self.a/2.)/3.), "lj")
        return vval

    def hcpcell0(self):
        return self.potential(self.distance(self.a/2., self.a/3./2.*(3.**0.5), 6.**0.5*self.a/3.))


mysys = System()
esc = []; efcc = []; ehcp = []
x = np.linspace(.95, 2., 50)
for d in x:
    mysys.a = d
    ef = 0.
    eh = 0.
    mysys.latticesc()
    esc.append(mysys.energy)
    mysys.latticefcc()
    efcc.append(mysys.energy)
    mysys.latticehcp()
    ehcp.append(mysys.energy)
    
# plt.plot(x, efcc, label='FCC')
# # plt.plot(x, ehcp, label='HCP')
# plt.grid()
# plt.legend()
# plt.xlabel('Lattice spacing [$\sigma$]')
# plt.ylabel('Energy [$\epsilon$]')
# plt.show()
print("sc min = ", min(esc), ":: fcc min = ", min(efcc), ":: hcp min = ", min(ehcp))