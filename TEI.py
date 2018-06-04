import math
import matplotlib.pyplot as plt


class var:
    # initialise static parameters
    dp1 = 7 * 100000  # Pa
    dp2 = 7 * 100000  # Pa
    mflow1 = 0.666 / 3
    mflow2 = 0.4799 / 3
    Rbarin1 = 0.8
    Rbarin2 = 3
    n1 = 2
    n2 = 6
    lbarn1 = 2
    lbarn2 = 2
    t_wall = 0.002  # m
    Km = mflow2 / mflow1
    alphaspray1 = 80  # DEG
    dr = 0.0003  # m
    # fluid properties LOX
    density1 = 1140.7  # kg/m^3
    kinvis1 = 1.736 * 10 ** (-7)  # Kg/ms

    # fluid properties Fuel
    density2 = 929.9  # kg/m^3
    kinvis2 = 2.122 * 10 ** (-7)  # Kg/ms

    # initialise dynamic parameters
    # stage 1
    alpha1 = 0
    phi1 = 0
    A1 = 0
    mu1 = 0
    Rn1 = 0
    Rin1 = 0
    Rein1 = 0
    rin1 = 0
    lin1 = 0
    ln1 = 0
    ls1 = 0
    R1 = 0
    rm = 0
    rm1 = 0
    factorlin1 = 4
    factorls1 = 2
    rmRnratio1 = 0

    # stage 2
    factorls2 = 2
    dp2g = [0]
    rm2 = []
    rin2 = []
    mu2 = []
    A2 = []
    Rn2 = []
    Rin2 = []
    Rein2 = 0
    alphaspray2 = 0
    phi2 = 0
    ln2 = 0
    ls2 = 0
    rmRnratio2 = []

    # combined
    spraytotal = 0
    tresidence = 0.001  # s
    lmix = 0
    deltalbarn = 0
    Rbarn = 0


class varint:
    # initialise static parameters
    dp1 = 7 * 100000  # Pa
    dp2 = 7 * 100000  # Pa
    mflow1 = 0.666 / 3
    mflow2 = 0.4799 / 3
    Rbarin1 = 0.8
    Rbarin2 = 3
    n1 = 2
    n2 = 3
    lbarn1 = 2
    lbarn2 = 2
    t_wall = 0.0008  # m
    Km = mflow2 / mflow1
    alphaspray1 = 70  # DEG
    # fluid properties LOX
    density1 = 1140.7  # kg/m^3
    kinvis1 = 1.736 * 10 ** (-7)  # Kg/ms

    # fluid properties Fuel
    density2 = 929.9  # kg/m^3
    kinvis2 = 2.122 * 10 ** (-7)  # Kg/ms

    # initialise dynamic parameters
    # stage 1
    alpha1 = 0
    phi1 = 0
    A1 = 0
    mu1 = 0
    Rn1 = 0
    Rin1 = 0
    Rein1 = 0
    rin1 = 0
    lin1 = 0
    ln1 = 0
    ls1 = 0
    R1 = 0
    rm = 0
    rm1 = 0
    factorlin1 = 4
    factorls1 = 2
    rmRnratio1 = 0

    # stage 2
    rm2 = []
    rin2 = 0
    mu2 = []
    A2 = []
    Rn2 = []
    Rin2 = 0
    Rein2 = 0
    alphaspray2 = 0
    phi2 = 0
    ln2 = 0
    rmRnratio2 = []

    # combined
    spraytotal = 0
    tresidence = 0.001  # s
    lmix = 0


def external1():
    var.alphaspray2 = (2 * var.alphaspray1 - 10) / 2
    var.A1 = (float(raw_input("A1. 2 alpha1 =     " + str(var.alphaspray1))))
    var.A2 = float(raw_input("A2. 2 alpha2 =     " + str(var.alphaspray2)))
    var.mu1 = float(raw_input("mu1. A1 =          " + str(var.A1)))
    var.mu2 = (float(raw_input("mu2. A2 =          " + str(var.A2))))

    # stage 1

    var.Rn1 = (0.475 * math.sqrt(var.mflow1 /
                                 (var.mu1 * math.sqrt(var.density1 * var.dp1))))
    var.Rin1 = (var.Rbarin1 * var.Rn1)
    var.rin1 = (math.sqrt(var.Rin1 * var.Rn1) / (var.n1 * var.A1))
    var.Rein1 = (2 * var.mflow1 /
                 (math.pi * math.sqrt(var.n1) * var.rin1 *
                  var.density1 * var.kinvis1))
    # stage 2
    var.Rn2 = (0.475 * math.sqrt(var.mflow2 /
                                 (var.mu2 * math.sqrt(var.density2 * var.dp2))))
    var.Rin2 = (var.Rbarin2 * var.Rn2)
    var.rin2 = (math.sqrt(var.Rin2 * var.Rn2) / (var.n2 * var.A2))
    var.Rein2 = (2 * var.mflow2 / (math.pi * math.sqrt(var.n2) * var.rin2 *
                                   var.density2 * var.kinvis2))

    if var.Rein1 < (10 ** 4):
        print "Re_in1 = " + str(var.Rein1) + "\nRe_in1 not big enough."
        return

    var.lin1 = var.factorlin1 * var.rin1
    var.ln1 = 2 * var.Rn1 * var.lbarn1
    var.ls1 = var.factorls1 * var.Rin1
    var.R1 = var.Rn1 + var.t_wall

    print("Re_in1 =  " + str(var.Rein1)
          + "\nA1 =      " + str(var.A1)
          + "\nRn1 =     " + str(var.Rn1)
          + "\nRin1 =    " + str(var.Rin1)
          + "\nrin1 =    " + str(var.rin1)
          + "\nln1 =    " + str(var.ln1)
          + "\nls1 =     " + str(var.ls1)
          + "\nR1 =      " + str(var.R1) + " (preliminary)"
          + "\nrm1 =     " + str(var.rm1)
          + "\n--== STAGE 1 COMPLETE ==--\n")

    if var.Rein2 < (10 ** 4):
        print "Re_in2 = " + str(var.Rein2) + "\nRe_in2 not big enough."
        return

    var.ln2 = 2 * var.Rn2 * var.lbarn2
    var.ls2 = var.factorls2 * var.Rin2

    print("Re_in2 =  " + str(var.Rein2)
          + "\nA2 =      " + str(var.A2)
          + "\nRn2 =     " + str(var.Rn2)
          + "\nRin2 =    " + str(var.Rin2)
          + "\nrin2 =    " + str(var.rin2)
          + "\nln2 =    " + str(var.ln2)
          + "\nls2 =     " + str(var.ls2)
          + "\n--== STAGE 2 COMPLETE ==--\n")


def external2():
    # combined
    var.alphaspray2 = (2 * var.alphaspray1 - 15) / 2
    var.A1 = (float(raw_input("A1. 2 alpha1 =     " + str(var.alphaspray1))))
    var.A2.append(float(raw_input("A2. 2 alpha2 =     " + str(var.alphaspray2))))
    var.mu1 = (float(raw_input("mu1. A1 =          " + str(var.A1))))
    var.mu2.append(float(raw_input("mu2. A2 =          " + str(var.A2[-1]))))

    # stage 1
    var.Rn1 = (0.475 * math.sqrt(var.mflow1 /
                                 (var.mu1 * math.sqrt(var.density1 * var.dp1))))
    var.Rin1 = (var.Rbarin1 * var.Rn1)
    var.rin1 = (math.sqrt(var.Rin1 * var.Rn1) / (var.n1 * var.A1))
    var.Rein1 = (2 * var.mflow1 /
                 (math.pi * math.sqrt(var.n1) * var.rin1 *
                  var.density1 * var.kinvis1))
    if var.Rein1 < (10 ** 4):
        print "Re_in1 = " + str(var.Rein1) + "\nRe_in1 not big enough."
        return

    var.lin1 = var.factorlin1 * var.rin1
    var.ln1 = 2 * var.Rn1
    var.ls1 = var.factorls1 * var.Rin1
    var.R1 = var.Rn1 + var.t_wall
    var.rm1 = var.Rn1

    print("Re_in1 =  " + str(var.Rein1)
          + "\nA1 =      " + str(var.A1)
          + "\nRn1 =     " + str(var.Rn1)
          + "\nRin1 =    " + str(var.Rin1)
          + "\nrin1 =    " + str(var.rin1)
          + "\nlin1 =    " + str(var.lin1)
          + "\nls1 =     " + str(var.ls1)
          + "\nR1 =      " + str(var.R1) + " (preliminary)"
          + "\nrm1 =     " + str(var.rm1)
          + "\n--== STAGE 1 COMPLETE ==--\n")

    # stage 2
    i = 0
    while abs(var.dp2g[i] - var.dp2) > 10000:
        print("Pressure drop difference = " + str(abs(var.dp2g[i] - var.dp2) / 100000) + " bar")
        if i >= 15:
            print("Nozzle too big.")
            return
        print var.R1
        var.Rn2.append(var.R1 + var.dr + (0.00025 * i))
        print var.Rn2[-1]
        var.deltalbarn = abs(1 / (var.lbarn2 / (2 * var.Rn2[-1])) * (var.lbarn1 /
                                                                     (2 * var.Rn1)))
        var.A2.append(float(raw_input("A2.  2 alpha2     =   "
                                      + str(var.alphaspray2)
                                      + "\n    deltal^bar_n2  =   "
                                      + str(var.deltalbarn))))

        var.Rbarn = abs(var.Rn2[-1] / var.Rn1)
        var.mu2.append(float(raw_input("mu2. R^bar_n2     =   "
                                       + str(var.Rbarn))))
        var.Rin2.append(var.Rbarin2 * var.Rn2[-1])
        print var.Rin2[-1]
        var.rin2.append(math.sqrt(var.Rin2[-1] * var.Rn2[-1] /
                                  (var.n2 * var.A2[-1])))
        print var.rin2[-1]
        var.dp2g.append(0.05 * var.mflow2 ** 2 /
                        (var.mu2[-1] ** 2 * var.density2 * var.Rn2[-1] ** 4))
        print var.dp2g

        i = i + 1


def stage1int():
    # det alpha, phi, A1, mu1
    """varint.A1 = float(raw_input("A. 2 alpha=" + varint.alphaspray1))
    varint.mu1=float(raw_input("Mu. A=" + varint.A1))
    varint.phi1=float(raw_input("Phi. A=" + varint.A1))
    varint.alpha=float((raw_input("Alpha. A="+ varint.A1))"""

    varint.A1 = 1.1
    varint.mu1 = 0.4
    varint.phi1 = 0.62
    varint.alpha = 50  # DEG

    varint.Rn1 = (0.475 * math.sqrt(varint.mflow1 /
                                    (varint.mu1 * math.sqrt(varint.density1 * varint.dp1))))
    varint.Rin1 = (varint.Rbarin1 * varint.Rn1)
    varint.rin1 = (math.sqrt(varint.Rin1 * varint.Rn1) / (varint.n1 * varint.A1))
    varint.Rein1 = (2 * varint.mflow1 /
                    (math.pi * math.sqrt(varint.n1) * varint.rin1 *
                     varint.density1 * varint.kinvis1))
    if varint.Rein1 < (10 ** 4):
        print "Re_in1 = " + str(varint.Rein1) + "\nRe_in1 not big enough."
        return

    varint.lin1 = varint.factorlin1 * varint.rin1
    varint.ln1 = 2 * varint.Rn1
    varint.ls1 = varint.factorls1 * varint.Rin1
    varint.R1 = varint.Rn1 + varint.t_wall
    varint.rm1 = varint.Rn1

    print("Re_in1 =  " + str(varint.Rein1)
          + "\nA1 =      " + str(varint.A1)
          + "\nRn1 =     " + str(varint.Rn1)
          + "\nRin1 =    " + str(varint.Rin1)
          + "\nrin1 =    " + str(varint.rin1)
          + "\nlin1 =    " + str(varint.lin1)
          + "\nls1 =     " + str(varint.ls1)
          + "\nR1 =      " + str(varint.R1)
          + "\nrm1 =     " + str(varint.rm1)
          + "\n--== STAGE 1 COMPLETE ==--\n")

    return


def stage2int():
    i = 0
    j = 0
    varint.rm2.append(varint.R1 + 0.0003)
    varint.Rn2.append(varint.rm2[i])
    varint.mu2.append(0.225 * varint.mflow2 / (varint.Rn2[i] ** 2 * math.sqrt(varint.density2 * varint.dp2)))

    # initialise iteration
    varint.A2.append(float(raw_input("A. mu_in=" + str(varint.mu2[i]))))
    varint.rmRnratio2.append(float(raw_input("r_m/R_m. A="
                                             + str(varint.A2[i]) + "& R^bar_in="
                                             + str(varint.Rbarin2))))
    varint.Rn2.append(varint.rm2[-1] / varint.rmRnratio2[-1])
    varint.mu2.append(0.225 * varint.mflow2 /
                      (varint.Rn2[-1] ** 2 * math.sqrt(varint.density2 * varint.dp2)))
    i = i + 2
    j = j + 1

    while (abs(varint.Rn2[-1] - varint.Rn2[-2]) > 0.0001):
        varint.A2.append(float(raw_input("A. mu_in=" + str(varint.mu2[-1]))))
        varint.rmRnratio2.append(float(raw_input("r_m/R_m. A="
                                                 + str(varint.A2[-1]) + "& R^bar_in="
                                                 + str(varint.Rbarin2))))
        varint.Rn2.append(varint.rm2[-1] / varint.rmRnratio2[-1])
        varint.mu2.append(0.225 * varint.mflow2 /
                          (varint.Rn2[-1] ** 2 * math.sqrt(varint.density2 * varint.dp2)))
        i = i + 2
        j = j + 1

    varint.Rin2 = varint.Rbarin2 * varint.Rn2[-1]
    varint.rin2 = math.sqrt(varint.Rin2 * varint.Rn2[-1] / (varint.n2 * varint.A2[-1]))
    varint.Rein2 = (2 * varint.mflow2 /
                    (math.pi * math.sqrt(varint.n2) * varint.rin2
                     * varint.density2 * varint.kinvis2))

    if varint.Rein2 < (10 ** 4):
        print "Re_in2 = " + str(varint.Rein2) + "\nRe_in1 not big enough."
        return

    if i > 15:
        print("Iteration of Rn2 does not converge.")
        return

    print("Re_in2 =  " + str(varint.Rein2)
          + "\nA2 =      " + str(varint.A2[-1])
          + "\nRn2 =     " + str(varint.Rn2[-1])
          + "\nRin2 =    " + str(varint.Rin2)
          + "\nrin2 =    " + str(varint.rin2)
          + "\nrm2 =     " + str(varint.rm2[-1])
          + "\n--== STAGE 2 COMPLETE ==--\n")
    return


def combinedint():
    varint.alphaspray2 = (float(raw_input("2 alpha. A2=" + str(varint.A2[-1]))))
    varint.spraytotal = varint.alphaspray2 - 35  # DEG
    varint.lmix = math.sqrt(2) * varint.tresidence * (varint.Km * varint.mu2[-1] /
                                                      ((varint.Km + 1) * (1 - varint.rmRnratio2[-1] ** 2))
                                                      * math.sqrt(varint.dp2 / varint.density2)
                                                      + varint.mu1 / (varint.Km + 1) * (1 - varint.rmRnratio1 ** 2)
                                                      * math.sqrt(varint.dp1 / varint.density1))
    varint.ln2 = 2 * varint.lbarn2 * varint.Rn2[-1]
    print("lmix =     " + str(varint.lmix)
          + "\nln2 =     " + str(varint.ln2)
          + "\ndeltaln = " + str(varint.ln2 - varint.lmix))


def resultsint():
    # stage 1
    s1 = [varint.Rin1 + varint.rin1,
          varint.Rin1,
          varint.Rn1,
          varint.R1]
    l1 = [varint.ls1 + varint.ln1,
          varint.ls1 + varint.ln1,
          0,
          0]
    plt.scatter(s1, l1)

    # stage 2
    s2 = [varint.Rin2 + varint.rin2,
          varint.Rin2,
          varint.Rn2[-1]]
    l2 = [varint.ln2,
          varint.ln2,
          0]
    plt.scatter(s2, l2, color="green")
    # combined
    plt.xlabel("Radius [m]")
    plt.ylabel("Height [m]")
    plt.title("Stage 1 and 2")
    plt.show()
    return


def resultsext():
    # stage 1
    s1 = [var.Rin1 + var.rin1,
          # var.Rin1,
          var.Rn1,
          var.Rn1,
          var.R1,
          var.R1,
          var.R1]
    l1 = [var.ls1 + var.ln1,
          # var.ls1+var.ln1,
          var.ln1,
          0,
          0,
          var.ln1,
          var.ls1 + var.ln1]

    # stage 2
    s2 = [var.Rin2 + var.rin2,
          # var.Rin2,
          var.Rn2,
          var.Rn2,
          0.01,
          0.01,
          0.01]
    l2 = [var.ls2 + var.ln2,
          # var.ls1+var.ln1,,
          var.ln2,
          0,
          0,
          var.ln2,
          var.ls2 + var.ln2]
    inletsr = [var.Rin1,
               var.Rin2]
    inletsh = [var.ls1 + var.ln1,
               var.ls2 + var.ln2]
    plt.scatter(inletsr, inletsh, color="red", label="inlet positions")
    plt.scatter(s1, l1, label="stage 1")
    plt.plot(s1, l1)
    plt.scatter(s2, l2, color="green", label="stage 2")
    plt.plot(s2, l2, color="green")
    plt.legend()

    # combined
    plt.xlabel("Radius [m]")
    plt.ylabel("Height [m]")
    plt.title("Injector elements: adapted geometries (preliminary)")
    plt.xlim(0, 0.015)
    plt.ylim(-0.001, 0.03)
    plt.show()
    return


def mainint():
    stage1int()
    stage2int()
    combinedint()
    resultsint()
    return


def mainext():
    external1()
    # external2()
    resultsext()
    return


mainext()
