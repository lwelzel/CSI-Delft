#<editor-fold desc="INIT">
import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, brentq
from scipy.optimize import OptimizeWarning
from scipy.interpolate import interp1d
from time import sleep
from sys import exit
# from pynverse import inversefunc
warnings.simplefilter("ignore", category=(RuntimeWarning, OptimizeWarning))
#</editor-fold>

#<editor-fold desc="GENERAL METHODS">
def results(s1,s2):
    #PRINT RESULTS IN CONSOLE
    for i in (s1,s2):
        for variable, value in vars(i).items():
            print(variable+":       "+value)
    """print("A1:         " + str(s1.A))
    print("A2:         " + str(s2.A))
    print("mu1:        " + str(s1.mu))
    print("mu2:        " + str(s2.mu))
    print("Rn1:        " + str(s1.Rn))
    print("Rn2:        " + str(s2.Rn))
    print("rin1:       " + str(s1.rin))
    print("rin2:       " + str(s2.rin))
    print("Re_in1:     " + str(s1.Re_in))
    print("Re_in2:     " + str(s2.Re_in))"""
    print("--=======================-")

    #PRINT RESULTS IN FIGURE

    #PRINT REGRESSION IN FIGURE

    #SAVE TO .txt OR QUIT
    res=input("Save results to .txt:    s\n"
              "Quit:                    q")
    if res==s:
        for i in (s1,s2):
            t=True
        #SAVE TO .txt
    else:
        exit()
    return
#</editor-fold>

#<editor-fold desc="GENERAL FUNCTIONS">
def inverse_root(x, a, b):
    return ((1 / ((1 + x) ** a)) + b)
def inv_inverse_root(x, a, b, y):
    return (y - ((1 / ((1 + x) ** a)) + b))
def cal_inv_inverse_root(y, para, datax):
    a = para[0]
    b = para[1]
    return brentq(inv_inverse_root, min(datax), 2 * max(datax), args=(a, b, y,), rtol=0.01)


def inverse_root_lin(x, a, b, c):
    return ((1 / ((1 + x) ** a)) + b * x + c)
def inv_inverse_root_lin(x, a, b, c, y):
    return (y - ((1 / ((1 + x) ** a)) + b * x + c))
def cal_inv_inverse_root_lin(y, para, datax):
    a = para[0]
    b = para[1]
    c = para[2]
    return brentq(inv_inverse_root_lin, min(datax), 2 * max(datax), args=(a, b, c, y,), rtol=0.01)


def log(x, a, b, c):
    return (a * (np.log(b * (x + c))))
def inv_log(x, a, b, c, y):
    return (y - (a * (np.log(b * (x + c)))))
def cal_inv_log(y, para, datax):
    a = para[0]
    b = para[1]
    c = para[2]
    return brentq(inv_log, min(datax), 1.1 * max(datax), args=(a, b, c, y,), rtol=0.01)
#</editor-fold>

#<editor-fold desc="DESIGN EQUATIONS">
def R_n(mf,mu,density,dp):
    return 0.475*np.sqrt(mf/(mu*np.sqrt(density*dp)))

def r_in(Rin,Rn,n,A):
    return np.sqrt(Rin*Rn/(n*A))

def Re_in(mf,n,rin,density,kinvis):
    return 0.637*mf/(np.sqrt(n)*rin*density*kinvis)

def interpol_A(alpha,lbar,para,data): #alpha, para=[Par.alpha2,Par.alpha5], data=[Data.Aalpha2,Data.Aalpha5]
    A1 = cal_inv_log(alpha, para[0], data[0])
    A2 = cal_inv_log(alpha, para[1], data[1])
    interpolate = interp1d([2., 0.5], [A1, A2], kind="linear", bounds_error=False,
                              fill_value='extrapolate')  # x = idicator, y = interp value
    return (interpolate(lbar))
#</editor-fold>

# DEFINE: Propellant stageing
def setup():  # SELECT: stage.oxidiser() or stage.fuel()
    s1 = Stage()
    s2 = Stage()
    s1.oxidiser()       #set stage propellant
    s2.fuel()           #set stage propellant

    # <editor-fold desc="SETUP">
    input = Parameters()
    input.Rbarin()

    s=(s1,s2)
    para=((input.alphas1,input.alphas2),
         (input.lbarn1,input.lbarn2),
         (input.dps1,input.dps2),
         (input.n1,input.n2),
         (input.Rbarin1,input.Rbarin2),
         ())
    for i in range(len(s)):
        s[i].alpha = para[0][i]
        s[i].lbar = para[1][i]
        s[i].dp = para[2][i]
        s[i].n = para[3][i]
        s[i].Rbarin = para[4][i]
    # </editor-fold>
    return (s1, s2, input)

#<editor-fold desc="DESIGN PROCESSES">
def external_v1(s1, s2, input):  # EXTERNAL MIXING, S1 IN GAS VORTEX OF S2
    # ALPHA S2
    s2.alpha = s1.alpha - input.spray_cone_angle_difference / 2

    a, c=True, 0
    while a and c<25:
        c+=1
        # A
        s1.A =interpol_A(s1.alpha,s1.lbar,[Par.alpha2,Par.alpha5],[Data.Aalpha2,Data.Aalpha5])
        s2.A = interpol_A(s2.alpha, s2.lbar, [Par.alpha2, Par.alpha5], [Data.Aalpha2, Data.Aalpha5])
        #print( "A1:         " + str(s1.A))
        #print( "A2:         " + str(s2.A))

        #MU
        s1.mu = inverse_root_lin(s1.A, *Par.mu)
        s2.mu = inverse_root_lin(s2.A, *Par.mu)
        #print( "mu1:        " + str(s1.mu))
        #print( "mu2:        " + str(s2.mu))

        #RN
        s1.Rn=R_n(s1.mf,s1.mu,s1.density,s1.dp)
        s2.Rn = R_n(s2.mf,s2.mu, s2.density, s2.dp)
        #print( "Rn1:        " + str(s1.Rn))
        #print( "Rn2:        " + str(s2.Rn))

        #rin
        s1.rin=r_in(s1.Rbarin,s1.Rn,s1.n,s1.A)
        s2.rin=r_in(s2.Rbarin,s2.Rn,s2.n,s2.A)
        #print( "rin1:       " + str(s1.rin))
        #print( "rin2:       " + str(s2.rin))

        #RE_IN
        s1.Re_in = Re_in(s1.mf,s1.n,s1.rin,s1.density,s1.kinvis)
        s2.Re_in = Re_in(s2.mf, s2.n, s2.rin, s2.density, s2.kinvis)
        #print( "Re_in1:     " + str(s1.Re_in))
        #print( "Re_in2:     " + str(s2.Re_in))

        if s1.Re_in > 10**4 and s2.Re_in > 10**4:
            a=False
    return s1,s2
def external_v2(s1, s2, input):  # EXTERNAL MIXING, S1 NOT IN GAS VORTEX OF S2
    #LRTC v2 pdf p.:74, also read at least 72 and 73
    t = True
    return s1,s2
def internal(s1, s2, input):  # INTERNAL MIXING
    #LRTC internal pdf p.:75-77, also read 72
    t = True
    return s1,s2
def complete(s1,s2,input):
    for i in (s1,s2):
        for variable, value in vars(i).items():
            if type(value)==str:
                t=True
                #working on atm
    return s1,s2


def main():
    s1, s2, input = setup()
    if input.external_mixing:
        s1, s2=internal(s1, s2, input)
    elif input.external_version1:
        s1,s2=external_v1(s1, s2, input)
    elif not (input.external_version1):
        s1, s2=external_v2(s1, s2, input)
    else:
        print( "Select proper design data")
    s1,s2=complete(s1,s2,input)
    results(s1,s2)
#</editor-fold>

# DEFINE: design, spray, geometry | if unknown put "na"
class Parameters:
    def __init__(self):
        # DESIGN
        self.external_mixing = False  # True or False
        self.external_version1 = True  # True for Version 1 False for Version 2
        self.dps1 = 7 * 10 ** 6  # Pa
        self.dps2 = 7 * 10 ** 6  # Pa

        # SPRAY
        self.spray_cone_angle_difference = 7.5  # DEG, between 10 and 15 DEG
        self.alphas1 = 80 /2  # DEG
        self.alphas2 = "na"

        # GEOMETRY
        self.Rbarin1_config = "open"  # open, closed
        self.Rbarin2_config = "closed"  # open, closed
        self.n1 = 3
        self.n2 = 6
        self.lbarn1 = 2.
        self.lbarn2 = 1.
        self.t_wall = 0.001  # m

        # GENERAL
        self.Re_cutoff = 10 ** 4

    def Rbarin(cls):
        if cls.Rbarin1_config == "closed":
            cls.Rbarin1 = 3.
        elif cls.Rbarin1_config == "open":
            cls.Rbarin1 = 0.75
        else:
            print( "Wrong Rbarin statement, put 'open' or 'closed'")

        if cls.Rbarin2_config == "closed":
            cls.Rbarin2 = 3.
        elif cls.Rbarin2_config == "open":
            cls.Rbarin2 = 0.75
        else:
            print( "Wrong Rbarin statement, put 'open' or 'closed'")

# DEFINE: Propellant properties
class Stage:
    def __init__(self, mf="NAN", density="NAN", kinvis="NAN",
                 dp="NAN",alpha="NAN",lbar="NAN",n="NAN",A="NAN",
                 mu="NAN",Rn="NAN",Rin="NAN", rin="NAN",Re_in="NAN",Rbarin="NAN"):
        self.mf = mf
        self.density = density
        self.kinvis = kinvis
        self.dp = dp
        self.alpha=alpha
        self.lbar=lbar
        self.n= n
        self.A = A
        self.mu = mu
        self.Rn = Rn
        self.Rin =Rin
        self.rin = rin
        self.Re_in = Re_in
        self.Rbarin= Rbarin


    def oxidiser(self):  # DEFINE: oxidiser properties
        self.mf = 0.666 / 3  # kg/s
        self.density = 1140.7  # kg/m^3
        self.kinvis = 1.736 * 10 ** (-7)  # Kg/ms
        return self

    def fuel(self):  # DEFINE: fuel properties
        self.mf = 0.4799 / 3  # kg/s
        self.density = 929.9  # kg/m^3
        self.kinvis = 2.122 * 10 ** (-7)  # Kg/ms
        return self

#<editor-fold desc="DATA">
class Data:
    mu = [1., .3, 0.9 * 0.2, 0.6 * 0.2, 0.16 * 0.2]
    Amu = [0., 2., 4., 6.4, 12]

    alpha2 = [0.5 * x for x in [60., 90., 100., 102.]]
    Aalpha2 = [1., 3., 7., 8.]

    alpha5 = [0.5 * x for x in [60., 92., 110., 113.]]
    Aalpha5 = [1., 2., 6., 8.]

    mu_in = [0.4, .25, .16, .12, 0.1]
    Amu_in = [.9, 2., 4., 6., 8.]

    rmRn1 = [0.3, 0.5, 0.69, 0.78]
    ArmRn1 = [.6, 1.4, 4., 8.]

    rmRn4 = [0.3, 0.4, 0.5, 0.6]
    ArmRn4 = [1.2, 2., 4., 8.]
class Par:
    mu, _ = curve_fit(inverse_root_lin, Data.Amu, Data.mu)
    alpha2, _ = curve_fit(log, Data.Aalpha2, Data.alpha2)
    alpha5, _ = curve_fit(log, Data.Aalpha5, Data.alpha5)
    muin, _ = curve_fit(inverse_root, Data.Amu_in, Data.mu_in)
    rmRn1, _ = curve_fit(log, Data.ArmRn1, Data.rmRn1)
    rmRn4, _ = curve_fit(log, Data.ArmRn4, Data.rmRn4)
#</editor-fold>

main()