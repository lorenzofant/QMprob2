import matplotlib.pyplot as plt
import numpy as np

class sistema():
    def __init__(self):
        self.a = 1.
        self.n = 16.
        self.epsilon = 1.
        self.sigma = 1.

    def potenziale(self, d):
        return 4.*self.epsilon*(pow(self.sigma/d, 12) - pow(self.sigma/d, 6))

    def distanza(self, x, y, z):
        d = (x**2 + y**2 + z**2)**0.5
        if d < self.a*1.1:
            self.closen += 1
        return d

    def latticesc(self):
        self.closen = 0.
        self.energia = self.sccell0()
        for i in np.linspace(-self.n, self.n, 2.*self.n+1):
            for j in np.linspace(-self.n, self.n, 2.*self.n+1):
                for k in np.linspace(-self.n, self.n, 2.*self.n+1):
                    if all(np.array([i, j, k]) == np.array([0., 0., 0.])):
                        c = 1
                    else:
                        self.energia += self.sccell(self.a*i, self.a*j, self.a*k)

    def latticefcc(self):
        self.closen = 0.
        #self.a = 1.5416
        self.energia = self.fcccell0()
        for i in np.linspace(-self.n, self.n, 2.*self.n+1):
            for j in np.linspace(-self.n, self.n, 2.*self.n+1):
                for k in np.linspace(-self.n, self.n, 2.*self.n+1):
                    if all(np.array([i, j, k]) == np.array([0., 0., 0.])):
                        c = 1
                    else:
                        self.energia += self.fcccell(self.a*i, self.a*j, self.a*k)

    def latticebcc(self):
        self.closen = 0.
        self.energia = self.bcccell0()
        for i in np.linspace(-self.n, self.n, 2.*self.n+1):
            for j in np.linspace(-self.n, self.n, 2.*self.n+1):
                for k in np.linspace(-self.n, self.n, 2.*self.n+1):
                    if all(np.array([i, j, k]) == np.array([0., 0., 0.])):
                        c = 1
                    else:
                        self.energia += self.bcccell(self.a*i, self.a*j, self.a*k)

    def latticehcp(self):
        self.closen = 0.
        #self.a = 1.09
        self.energia = self.hcpcell0()
        for i in np.linspace(-self.n, self.n, 2.*self.n+1.):
            for j in np.linspace(-self.n, self.n, 2.*self.n+1.):
                for k in np.linspace(-self.n/2., self.n/2., self.n+1.):
                    if all(np.array([i, j, k]) == np.array([0., 0., 0.])):
                        c = 1
                    else:
                        self.energia += self.hcpcell(self.a*i, self.a*j, self.a*k)

    def sccell(self, x, y, z):
        return self.potenziale(self.distanza(x, y, z))

    def fcccell(self, x, y, z):
        a = 2.**0.5
        return self.potenziale(self.distanza(a*x+(y+z)/a + self.noise*np.random.normal(), y/a + self.noise*np.random.normal(), z/a + self.noise*np.random.normal()))

    def bcccell(self, x, y, z):
        a = 1./(3.**0.5)
        return self.potenziale(self.distanza(2.*a*x, 2.*a*y, 2.*a*z)) + self.potenziale(self.distanza(2.*a*x+self.a*a, 2.*a*y+self.a*a, 2.*a*z+self.a*a))

    def hcpcell(self, x, y, z):
        return self.potenziale(self.distanza(x + y/2. + self.noise*np.random.normal(), y/2.*(3.**0.5)+self.noise*np.random.normal(), 24.**0.5*z/3.+self.noise*np.random.normal())) + self.potenziale(self.distanza(x+(y+self.a)/2.+self.noise*np.random.normal(), (y+self.a/3.)/2.*(3.**0.5)+self.noise*np.random.normal(), 24.**0.5*(z+self.a/2.)/3.+self.noise*np.random.normal()))


    def sccell0(self):
        return 0.

    def fcccell0(self):
        return 0.


    def bcccell0(self):
        a = 1./2.**0.5
        return self.potenziale(self.distanza(self.a*a, self.a*a, self.a*a))

    def hcpcell0(self):
        return self.potenziale(self.distanza(self.a/2.+self.noise*np.random.normal(), self.a/3./2.*(3.**0.5)+self.noise*np.random.normal(), 6.**0.5*self.a/3.+self.noise*np.random.normal()))


a = sistema()
esc = []
efcc = []
ebcc = []
ehcp = []
x = np.linspace(.95, 1.2, 41)
for d in x:
    print d
    a.noise = 0.
    a.a = d
    ef = 0.
    eh = 0.
    a.latticesc()
    esc.append(a.energia)
    for i in range(1):
        #print i
        #a.a = 1.54179
        a.latticefcc()
        ef += (a.energia)
        #a.a = 1.0902
        a.latticehcp()
        eh += (a.energia)
    efcc.append(ef/1.)
    print efcc[-1]
    a.latticebcc()
    ebcc.append(a.energia)
    ehcp.append(eh/1.)
    print ehcp[-1]
plt.plot(x, esc, label='SC')
plt.plot(x, efcc, label='FCC')
print sum(efcc)/len(efcc)
plt.plot(x, ebcc, label='BCC')
print sum(ehcp)/len(ehcp)
plt.plot(x, ehcp, label='HCP')
plt.grid()
plt.legend()
plt.xlabel('Lattice spacing [$\sigma$]')
plt.ylabel('Energy [$\epsilon$]')
plt.show()
np.savetxt('energie7', [x, efcc, ehcp])
