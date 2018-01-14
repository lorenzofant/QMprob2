import matplotlib.pyplot as plt
import numpy as np
import itertools

class System():
    def __init__(self):
        self.a = 1.0902#3.75
        self.n = 4.
        # constants for LJ potential
        self.epsilon = 1. #137.37 #99.55
        self.sigma = 1. #3.4
        self.mass = 1.
        # constants for BFW potential
        self.bfw_epsilon=142.1; self.bfw_sigma=3.76
        self.bfw_C6 = -1.10727; self.bfw_A0 = 0.27783; self.bfw_A3 = -25.2696
        self.bfw_C8 = -0.16971; self.bfw_A1 = -4.50431; self.bfw_A4 = -102.0195
        self.bfw_alpha = 12.5; self.bfw_C10 = -0.01361; self.bfw_A2 = -8.33122; self.bfw_A5 = -113.25
        self.bfw_delta = 0.01
        # constants for BBMS potential
        self.bbms_epsilon=142.1; self.bbms_sigma=3.76
        self.bbms_C6 = -1.10727; self.bbms_A0 = 0.27783; self.bbms_A3 = -25.2696
        self.bbms_C8 = -0.16971; self.bbms_A1 = -4.50431; self.bbms_A4 = -102.0195
        self.bbms_alpha = 12.5; self.bbms_C10 = -0.01361; self.bbms_A2 = -8.33122; self.bbms_A5 = -113.25
        self.bbms_delta = 0.01
        # constants for hdf potential
        self.hfd_epsilon = 143.224;self.hfd_C6 = -1.0914254;self.hfd_A = 0.9502720e7
        self.hfd_sigma =3.759; self.hfd_C8 = -0.6002595
        self.hfd_alpha = 16.345655; self.hfd_C10 = -0.3700113; self.hfd_gamma = 2.0

    def potential(self, d, vtype = "lj"):
        vval=0.
        if vtype=="lj": 
            r = d/self.sigma # to avoid recalculation of this term when calling pow()
            vval= 4*self.epsilon*(pow(r, -12) - pow(r, -6))

        if vtype=="bfw": 
            r = d/self.bfw_sigma
            vval= self.bfw_epsilon*(np.exp(self.bfw_alpha*(1-r))*(self.bfw_A1*(r-1)+self.bfw_A2*pow((r-1),2)+self.bfw_A3*pow((r-1),3)+self.bfw_A4*pow((r-1),4)+self.bfw_A5*pow((r-1),5))+(self.bfw_C6/(self.bfw_delta+pow(r,6))+self.bfw_C8/(self.bfw_delta+pow(r,8))+self.bfw_C10/(self.bfw_delta+pow(r,10))))
            
        if vtype=="bbms": 
            r = d/self.bbms_sigma
            vval= self.bbms_epsilon*(np.exp(self.bbms_alpha*(1-r))*(self.bbms_A1*(r-1)+self.bbms_A2*pow((r-1),2)+self.bbms_A3*pow((r-1),3)+self.bbms_A4*pow((r-1),4)+self.bbms_A5*pow((r-1),5))+(self.bbms_C6/(self.bbms_delta+pow(r,6))+self.bbms_C8/(self.bbms_delta+pow(r,8))+self.bbms_C10/(self.bbms_delta+pow(r,10)))+self.bbms_alpha*np.exp(-50.0*pow(r-1.33,2)))

        if vtype=="hfd": 
            r = d/self.hfd_sigma
            if r <= 1.4:
                vval = self.hfd_epsilon*(self.hfd_A*np.exp(-1*self.hfd_alpha*r)+(self.hfd_C6*pow(r,-6)+self.hfd_C8*pow(r,-8)+self.hfd_C10*pow(r,-10))*np.exp(-1*(pow(1.4/r,2))))
            else:
                vval = self.hfd_epsilon*(self.hfd_A*np.exp(-1*self.hfd_alpha*r)+(self.hfd_C6*pow(r,-6)+self.hfd_C8*pow(r,-8)+self.hfd_C10*pow(r,-10)))
        return vval

    def potentialp(self, x, i, j):
        d = np.linalg.norm(x)
        return self.LJpotential_deriv2(d)*np.dot(x,i)*np.dot(x,j)/(d**2) + self.LJpotential_deriv(d)*(np.dot(i,j)/d - np.dot(x,i)*np.dot(x,j)/(d**3))

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
            c = np.array([0,0,0])
            c[l]+=1
            a = []
            for m in range(3):
                d = np.array([0.,0.,0.])
                d[m]+=1.
                self.energy = 0.#self.hcpcell0(potentialType)
                for i, j, k in itertools.product(r,r,r):
                    x = np.array([self.a*i, self.a*j, self.a*k])
                    if (i,j,k) != (0,0,0): self.energy += self.fcccellp(x, c, d)*self.phasek(x)
                self.energy /= self.mass
                a.append(self.energy)
            self.dm.append(a)
        self.w = np.linalg.eigvals(self.dm)

    def phasek(self, x):
        return np.e**(1.j*np.dot(self.k,x))

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


    def hcpcellp(self, x, c, d):
        return self.potentialp([x[0] + x[1]/2., x[1]/2.*(3.**0.5), 24.**0.5*x[2]/3.], c, d) + self.potentialp([x[0] + (x[1] + self.a)/2., (x[1] + self.a/3.)/2.*(3.**0.5), 24.**0.5*(x[2] + self.a/2.)/3.], c, d)

    def fcccellp(self, x, c, d):
        a = 2.**0.5
        return self.potentialp([a*x[0]+(x[1]+x[2])/a, x[1]/a, x[2]/a], c, d)


# main body of the program
a = System()
esc = []; ebcc = []
efcc_lj = []; efcc_bfw = []; efcc_hfd = []; efcc_bbms = []
ehcp_lj = []; ehcp_bfw = []; ehcp_hfd = []; ehcp_bbms = []
x = np.linspace(0.,(3./2.)**0.5*np.pi,100)
w = []
for el in x:
    print el
    a.k = el*np.array([0.,0.,1.])
    a.dynamicmatrix()
    w.append(a.w)
w = np.transpose(w)
for el in w:
    plt.plot(x,el)
plt.show()
x = np.linspace(1.03,1.12,20)
a.n =14.
for d in x:
    print( d )
    a.a = d
    ef = 0.
    eh = 0.
    a.lattice("fcc", "lj")
    efcc_lj.append(a.energy)
    a.lattice("hcp", "lj")
    ehcp_lj.append(a.energy)

    # a.lattice("fcc", "bfw")
    # efcc_bfw.append(a.energy)
    # a.lattice("hcp", "bfw")
    # ehcp_bfw.append(a.energy)

    # a.lattice("fcc", "hfd")
    # efcc_hfd.append(a.energy)
    # a.lattice("hcp", "hfd")
    # ehcp_hfd.append(a.energy)

    # a.lattice("fcc", "bbms")
    # efcc_bbms.append(a.energy)
    # a.lattice("hcp", "bbms")
    # ehcp_bbms.append(a.energy)

    # a.lattice("sc", "hfd")
    # esc.append(a.energy)
    # for i in range(1):
    #     #print( i )
    #     #a.a = 1.54179
    #     a.lattice("fcc", "bfw")
    #     ef += (a.energy)
    #     #a.a = 1.0902
    #     a.lattice("hcp", "bfw")
    #     eh += (a.energy)
    # efcc.append(ef/1.)
    # print( efcc[-1] )
    # a.lattice("bcc", "hfd")
    # ebcc.append(a.energy)
    # ehcp.append(eh/1.)
    # print( ehcp[-1] )
# plt.plot(x, esc, label='SC')

plt.plot(x, efcc_lj, label='FCC with LJ')
# plt.plot(x, efcc_bfw, label='FCC with BFW')
# plt.plot(x, efcc_hfd, label='FCC with HFD')
# plt.plot(x, efcc_bbms, label='FCC with BBMS')

# print( sum(efcc)/len(efcc) )
# plt.plot(x, ebcc, label='BCC')
# print( sum(ehcp)/len(ehcp) )
plt.plot(x, ehcp_lj, label='HCP with LJ')
# plt.plot(x, ehcp_bfw, label='HCP with BFW')
# plt.plot(x, ehcp_hfd, label='HCP with HFD')
# plt.plot(x, ehcp_bbms, label='HCP with BBMS')
plt.grid()
plt.legend()
plt.xlabel('Lattice spacing [$\sigma$]')
plt.ylabel('Energy [$\epsilon$]')
plt.show()
# np.savetxt('energie6', [x, efcc, ehcp])
# print("LJ, fcc min = ", min(efcc_lj), ":: hcp min = ", min(ehcp_lj), ":: difference factor ",1-(min(efcc_lj)/min(ehcp_lj)))
# print("BFW, fcc min = ", min(efcc_bfw), ":: hcp min = ", min(ehcp_bfw), ":: difference factor ",1-(min(efcc_bfw)/min(ehcp_bfw)))
# print("BFW, fcc min = ", min(efcc_hfd), ":: hcp min = ", min(ehcp_hfd), ":: difference factor ",1-(min(efcc_hfd)/min(ehcp_hfd)))
# print("BFW, fcc min = ", min(efcc_bbms), ":: hcp min = ", min(ehcp_bbms), ":: difference factor ",1-(min(efcc_bbms)/min(ehcp_bbms)))
