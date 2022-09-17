import math
import matplotlib.pyplot as plt
import numpy as np


def calculatePh(H_concentration):
    return -math.log10(H_concentration)


def calculateH_concentration(Ph):
    return math.pow(10, -Ph)


def convertPh_to_POH(ph_POH):
    return 14 - ph_POH


def convertKa_to_pka(ka):
    return -math.log10(ka)


def convertPka_to_ka(pka):
    return math.pow(10, -pka)


def convertPka_to_Pkb(Pka):
    return 14 - Pka


def convertPka_to_ka(pka):
    return math.pow(10, -pka)


def convertKa_to_Kb(Ka1):
    AUTOIONIZATION_OF_WATER = 1e-14
    return AUTOIONIZATION_OF_WATER / Ka1


def normality(molarity, mole_ratio):
    normality = molarity * mole_ratio
    return normality


def quadSolver_solutionVariant(ka, c):
    discriminant = math.sqrt(ka**2 + (4 * ka * c))
    ph = -math.log10((discriminant - ka) / 2)
    return ph

# henderson-Hasselbalch equations


def Henderson_Hasselbalch(Pka, A_minus, HA):
    return Pka + math.log10(A_minus / HA)


def allam(ka, A_minus, HA):
    return math.log10((A_minus / HA) / ka)

# extra equations


def acidBase_equilibriumCalculator():
    print("Enter 1 for first disassociation ")
    print("Enter 2 for second disassociation ")
    disassociation_decider = int(input())
    if disassociation_decider == 1:
        initial_concentration = float(
            input("\nEnter initial Concentration of the solution: "))
        ka = float(input("Enter Ka of the solution: "))
        H_concentration = quadSolver_solutionVariant(ka, initial_concentration)
        print("The H+ Concentration of this solution is ", H_concentration, "\nThe ph is ", calculatePh(H_concentration))

    elif disassociation_decider == 2:
        A = 1
        initial_concentration = float(
            input("\nEnter initial Concentration of the solution: "))
        ka1 = float(input("Enter Ka1 of the solution: "))
        ka2 = float(input("Enter Ka2 of the solution: "))
        H_initialConcentration = quadSolver(A, ka1,(-ka1 * initial_concentration))
        H_2ndConcentration = quadSolver(A, (ka2 + H_initialConcentration),(-ka2 * H_initialConcentration))
        H_combinedConcentration = H_initialConcentration + H_2ndConcentration
        print("\nThe combined H+ concentration for this solution is ",  H_combinedConcentration, "\nThe combined ph is ",  calculatePh(H_combinedConcentration),
              "\nThe A- ion concentration is ", H_2ndConcentration)


def quadSolver(a, b, c):
    discriminant = (b**2) - (4 * a * c)
    x1 = (-b + math.sqrt(discriminant)) / (2 * a)
    x2 = (-b - math.sqrt(discriminant)) / (2 * a)
    if (x1 > 0):
        return x1
    else:
        return x2

# Functions used by titrations


def titrations():
    print("Enter one of the options below: ")
    print("\n1. Strong acid - Strong base: ")
    print("2. Strong base - Strong acid: ")
    print("3. Weak acid - Strong base: ")
    print("4. Weak base - Strong acid: ")
    titrations_decider = int(input())
    if (titrations_decider == 1):
        info1 = generate_prompt(titrations_decider)
        acid_mole_ratio = float(
            input("\nEnter number of H+ the acid can give up: "))
        base_mole_ratio = float(
            input("Enter number of H+ the base can receive: "))
        moles_of_acid = (info1[0] * acid_mole_ratio) * (info1[2] / 1000)
        numbers_of_volumes = parser(info1[3])
        list_of_phs = []
        midpoint_ph = 0
        for i in numbers_of_volumes:
            moles_of_base = (info1[1] * base_mole_ratio) * (i / 1000)
            total_volume = (i + info1[2]) / 1000
            check = moles_of_acid - moles_of_base
            # previous to eqpoint
            if check > 0:
                if (check == moles_of_acid):
                    conc_sol = (moles_of_acid - moles_of_base) / total_volume
                    midpoint_ph = -math.log10(conc_sol)
                    list_of_phs.append(midpoint_ph)
                conc_sol = (moles_of_acid - moles_of_base) / total_volume
            # Eq point
            elif (check == 0):
                conc_sol = 1.00e-7
            # Excess
            else:
                oh = (moles_of_base - moles_of_acid) / total_volume
                conc_sol = 1.00e-14 / oh

            list_of_phs.append(-math.log10(conc_sol))
        for i in range(len(list_of_phs)):
            print("\n", numbers_of_volumes[i], "ml. ", list_of_phs[i])
        print("\nMidpoint is at ph of ", midpoint_ph)
        print("Equivalency point is at ph of ", 7)
        Graph = int(input("\nEnter 1 to graph the titration curve: "))
        if Graph == 1:
            graph(numbers_of_volumes, list_of_phs, titrations_decider)

    elif (titrations_decider == 2):
        info1 = generate_prompt(titrations_decider)
        acid_mole_ratio = float(
            input("\nEnter number of H+ the acid can give up: "))
        base_mole_ratio = float(
            input("Enter number of H+ the base can receive: "))
        moles_of_base = (info1[0] * base_mole_ratio) * (info1[2] / 1000)
        numbers_of_volumes = parser(info1[3])
        list_of_phs = []
        midpoint_ph = 0
        for i in numbers_of_volumes:
            moles_of_acid = (info1[1] * acid_mole_ratio) * (i / 1000)
            total_volume = (i + info1[2]) / 1000
            check = moles_of_base - moles_of_acid
            # previous to eqpoint
            if check > 0:
                conc_sol = (moles_of_base - moles_of_acid) / total_volume
                if (check == moles_of_base):
                    conc_sol = (moles_of_base - moles_of_acid) / total_volume
                    midpoint_ph = 14 + math.log10(conc_sol)
                    list_of_phs.append(midpoint_ph)
                list_of_phs.append(14 + math.log10(conc_sol))
            # Eq point
            elif (check == 0):
                list_of_phs.append(7)
            # Excess
            else:
                h = (moles_of_acid - moles_of_base) / total_volume
                list_of_phs.append(-math.log(h))

        for i in range(len(list_of_phs)):
            print("\n", numbers_of_volumes[i], "ml. ", list_of_phs[i])
        print("\nMidpoint is at ph of ", midpoint_ph)
        print("Equivalency point is at ph of ", 7)
        Graph = int(input("\nEnter 1 to graph the titration curve: "))
        if Graph == 1:
            graph(numbers_of_volumes, list_of_phs, titrations_decider)

    elif (titrations_decider == 3):
        info2 = generate_prompt(titrations_decider)
        number_of_volumes = parser(info2[4])
        list_of_phs = []
        moles_of_acid = info2[1] * (info2[3] / 1000)
        midpoint_ph = 0
        eqpoint_ph = 0
        for i in number_of_volumes:
            total_volume = (i+info2[3])/1000
            moles_of_base = info2[1] * (i / 1000)
            check = moles_of_acid - moles_of_base
            # Starting off ml
            if (i == 0):
                ph = quadSolver_solutionVariant(info2[2], info2[0])
                list_of_phs.append(ph)
            # Buffer zone
            elif (check > 0):
                reaction_ratioCoefficient = (moles_of_acid - moles_of_base) / (moles_of_base)
                if (check == moles_of_acid):
                    midpoint_ph = 14 + \
math.log10(info2[2] * reaction_ratioCoefficient)
                    list_of_phs.append(midpoint_ph)
                ph = -math.log10(info2[2] * reaction_ratioCoefficient)
                list_of_phs.append(ph)
            # eq-point
            elif (check == 0):
                kb = convertKa_to_Kb(info2[2])
                C = info2[0]
                poh = quadSolver_solutionVariant(kb, C)
                eqpoint_ph = 14-poh
                list_of_phs.append(eqpoint_ph)
            # After eq-point excess
            elif (check < 0):
                ratio = ((moles_of_base) - (moles_of_acid)) / (total_volume)
                ph = 14 + math.log10(ratio)
                list_of_phs.append(ph)
        for i in range(len(list_of_phs)):
            print("\n", number_of_volumes[i], "ml. ", list_of_phs[i])
            
        print("\nMidpoint is at ph of ", midpoint_ph)
        print("Equivalency point is at ph of ", eqpoint_ph)
        Graph = int(input("\nEnter 1 to graph the titration curve: "))
        if Graph == 1:
            graph(number_of_volumes, list_of_phs, titrations_decider)

    elif titrations_decider == 4:
        info2 = generate_prompt(titrations_decider)
        number_of_volumes = parser(info2[4])
        list_of_phs = []
        moles_of_base = info2[1] * (info2[3] / 1000)
        midpoint_ph = 0
        eqpoint_ph = 0
        for i in number_of_volumes:
            total_volume = (i+info2[3])/1000
            moles_of_acid = info2[1] * (i / 1000)
            check = moles_of_base - moles_of_acid
            # Starting off ml
            if (i == 0):
                ph = 14-quadSolver_solutionVariant(info2[2], info2[0])
                list_of_phs.append(ph)
            # Buffer zone
            elif (check > 0):
                reaction_ratioCoefficient = (moles_of_base - moles_of_acid) / (moles_of_acid)
                if (check == moles_of_base):
                    midpoint_ph = 14+math.log10(info2[2] * reaction_ratioCoefficient)
                    list_of_phs.append(midpoint_ph)
                ph = 14+math.log10(info2[2] * reaction_ratioCoefficient)
                list_of_phs.append(ph)
            # eq-point
            elif (check == 0):
                C = info2[0]
                ka = convertKa_to_Kb(info2[2])
                eqpoint_ph = quadSolver_solutionVariant(ka, C)
                list_of_phs.append(eqpoint_ph)
            # After eq-point excess
            elif (check < 0):
                ratio = ((moles_of_acid) - (moles_of_base))/total_volume
                ph = -math.log10(ratio)
                list_of_phs.append(ph)
        for i in range(len(list_of_phs)):
            print("\n", number_of_volumes[i], "ml. ", list_of_phs[i])
            
        print("\nMidpoint is at ph of ", midpoint_ph)
        print("Equivalency point is at ph of ", eqpoint_ph)
        Graph = int(input("\nEnter 1 to graph the titration curve: "))
        if Graph == 1:
            graph(number_of_volumes, list_of_phs, titrations_decider)


def generate_prompt(decider):
    if decider == 1:
        conc_acid = float(input("\nEnter molarity of Strong acid: "))
        conc_base = float(input("Enter molarity of Strong base: "))
        volume_acid = float(input("\nEnter milliliters of Strong acid: "))
        list_Vb = str(input("Enter volume of Strong base ex:(1.0,2.0) "))
        return [conc_acid, conc_base, volume_acid, list_Vb]
    elif decider == 2:
        conc_base = float(input("\nEnter molarity of Strong base: "))
        conc_acid = float(input("Enter molarity of Strong acid: "))
        volume_base = float(input("\nEnter milliliters of Strong base: "))
        list_Va = str(input("Enter volume of Strong acid ex:(1.0,2.0) "))
        return [conc_base, conc_acid, volume_base, list_Va]

    elif decider == 3:
        conc_acid = float(input("\nEnter molarity of Weak acid : "))
        conc_base = float(input("Enter molarity of Strong base : "))
        ka = float(input("\nEnter ka : "))
        volume_acid = float(input("\nEnter milliliters of Weak acid: "))
        list_Vb = str(input("Enter volume of Strong base ex:(1.0,2.0) "))
        return [conc_acid, conc_base, ka, volume_acid, list_Vb]
    elif decider == 4:
        conc_acid = float(input("\nEnter molarity of Weak base : "))
        conc_base = float(input("Enter molarity of Strong acid : "))
        kb = float(input("\nEnter kb : "))
        volume_base = float(input("\nEnter milliliters of Weak base: "))
        list_Va = str(input("Enter volume of Acid base ex:(1.0,2.0) "))
        return [conc_base, conc_acid, kb, volume_base, list_Va]


def parser(list_Vb):
    list_of_volumes = list_Vb.split(',')
    numbers_of_volumes = []
    for i in list_of_volumes:
        volume = float(i)
        numbers_of_volumes.append(volume)
    return numbers_of_volumes


def graph(list_of_volumes, list_of_phs, decider):

    if decider == 1:
        xvalues = np.asarray(list_of_volumes)
        yvalues = np.asarray(list_of_phs)
        plt.plot(xvalues, yvalues)
        plt.title("Titration of Strong acid with Strong base")
        plt.xlabel("Volume of Base added in ml")
        plt.ylabel("Ph")
        plt.grid()
        plt.show()
    elif decider == 2:
        xvalues = np.asarray(list_of_volumes)
        yvalues = np.asarray(list_of_phs)
        plt.plot(xvalues, yvalues)
        plt.title("Titration of Strong base with Strong acid")
        plt.xlabel("Volume of Acid added in ml")
        plt.ylabel("Ph")
        plt.grid()
        plt.show()

    elif decider == 3:
        xvalues = np.asarray(list_of_volumes)
        yvalues = np.asarray(list_of_phs)
        plt.plot(xvalues, yvalues)
        plt.title("Titration of Weak acid with Strong base")
        plt.xlabel("Volume of Base added in ml")
        plt.ylabel("Ph")
        plt.grid()
        plt.show()
    elif decider == 4:
        xvalues = np.asarray(list_of_volumes)
        yvalues = np.asarray(list_of_phs)
        plt.plot(xvalues, yvalues)
        plt.title("Titration of Weak base with Strong acid")
        plt.xlabel("Volume of Acid added in ml")
        plt.ylabel("Ph")
        plt.grid()
        plt.show()

# piece Functions


def GeneralEquilibrium():

    print("\nGeneral Equilibrium beta 1.0v")
    print("Select one of the equations below")
    print("\n1. 1a + 1b --> 2c")
    print("2. 1a + 1b --> 1c + 1d")
    equation_decider = int(input())
    if(equation_decider == 1):
        Keq = float(input("\nEnter the Keq value "))
        a = float(input("Enter the Concentration of a "))
        b = float(input("Enter the Concentration of b "))
        if(Keq > 100):
            print("The reaction is a Large K reaction ")
            x = quadSolver(4-Keq, (-4*b + Keq*b - Keq*a), 4*b**2)
            print(x)
            # Finish the Large K reaction later
        elif(Keq > 0.01):
            print("\nThe reaction is Intermediate K reaction ")
            x = quadSolver(Keq-4, -Keq*(a+b), Keq*(a*b))
            equilibriumConcentration_a = a-x
            equilibriumConcentration_b = b-x
            equilibriumConcentration_c = 2*x
            print("\nThe x is: ", x)
            print("The Equilibrium Concentration of a is: ",
                  equilibriumConcentration_a)
            print("The Equilibrium Concentration of b is: ",
                  equilibriumConcentration_b)
            print("The Equilibrium Concentration of c is: ",
                  equilibriumConcentration_c)
        elif(Keq < 0.01):
            print("The reaction is Small K reaction ")
            x = quadSolver(Keq-4, -Keq*(a+b), Keq*(a*b))
            equilibriumConcentration_a = a-x
            equilibriumConcentration_b = b-x
            equilibriumConcentration_c = 2*x
            print("\nThe x is: ", x)
            print("The Equilibrium Concentration of a is: ",
                  equilibriumConcentration_a)
            print("The Equilibrium Concentration of b is: ",
                  equilibriumConcentration_b)
            print("The Equilibrium Concentration of c is: ",
                  equilibriumConcentration_c)


def Solution_main():
    check = True
    while (check == True):
        reader = SolutionMenu()
        if reader == 1:
            ka = float(input("Enter the ka value: "))
            Pka = convertKa_to_pka(ka)
            print("\nThe Pka is ", Pka, "\n")
        elif reader == 2:
            Pka = float(input("Enter the Pka "))
            ka = convertPka_to_ka(Pka)
            print("\nThe ka value is ", ka, "\n")
        elif reader == 3:
            Pka_Pkb = float(input("Enter the Pka or Pkb "))
            Pkb_Pka = convertPka_to_Pkb(Pka_Pkb)
            print("\nThe Poh or ph is", Pkb_Pka, "\n")
        elif reader == 4:
            pass
        elif reader == 5:
            titrations()
        elif reader == 6:
            print("\nEnding the program in 3..2..1 ")
            check = False


def SolutionMenu():

    print(
        "**************************************\n\tSolution PROGRAM\n\tBy: Vincent Allam\n**************************************"
    )
    "\n"
    print("Basic operations: ")
    print("\tOption 1: Convert ka to Pka ")
    print("\tOption 2: Convert Pka to Ka ")
    print("\tOption 3: Convert Pka to Pkb or vice versa ")
    print("\nAdvanced operations: ")
    print("\tOption 4: Common ion / buffer Problems ")
    print("\tOption 5: Titration Calculator ")
    print("Option 6: Exit Program")
    option_decider = int(input())
    return option_decider


def acid_base_main():
    check = True
    while (check == True):
        reader = acid_base_Menu()
        if reader == 1:
            H_concentration = float(input("Enter the H+ Concentration: "))
            ph = calculatePh(H_concentration)
            print("\nThe PH is ", ph, "\n")
        elif reader == 2:
            ph = float(input("Enter the Ph "))
            H_concentration = calculateH_concentration(ph)
            print("\nThe H+ Concentration is ", H_concentration, "\n")
        elif reader == 3:
            ph_PoH = float(input("Enter the ph or POH "))
            PoH_ph = convertPh_to_POH(ph_PoH)
            print("\nThe Poh or ph is", PoH_ph, "\n")
        elif reader == 4:
            acidBase_equilibriumCalculator()
        elif reader == 5:
            ka = float(input("Enter the ka value "))
            kb = convertKa_to_Kb(ka)
            print("\nThe conjugate kb value is ", kb)
        elif reader == 6:
            ka = float(input("Enter the ka value "))
            pka = convertKa_to_pka(ka)
            print("\nThe Pka is ", pka, "\n")
        elif reader == 7:
            pka = float(input("Enter the Pka value "))
            ka = convertPka_to_ka(pka)
            print("\nThe ka is ", ka, "\n")
        elif reader == 8:
            print("\nEnding the program in 3..2..1 ")
            check = False


def acid_base_Menu():

    print(
        "**************************************\n\tACID - BASE PROGRAM\n\tBy: Vincent Allam\n**************************************"
    )
    "\n"
    print("Basic operations: ")
    print("\tOption 1: Convert to Ph ")
    print("\tOption 2: Convert Ph to H+ concentration ")
    print("\tOption 3: Convert Ph to POh or vice versa ")
    print("\nAdvanced operations: ")
    print("\tOption 4: Acid-Base Equilibrium calculator")
    print("\tOption 5: Convert Ka to conjugate Kb ")
    print("\tOption 6: Convert Ka/Kb to PKa/Pkb ")
    print("\tOption 7: Convert Pka/Pkb to Ka/Kb ")
    print("Option 8: Exit Program")
    option_decider = int(input())
    return option_decider


def redox(any):
    return "red"


def generate_prompt_txt(decider):
    file1 = open()
    if decider == 1:
        conc_acid = float(input("\nEnter molarity of Strong acid: "))
        conc_base = float(input("Enter molarity of Strong base: "))
        volume_acid = float(input("\nEnter milliliters of Strong acid: "))
        list_Vb = str(input("Enter volume of Strong base ex:(1.0,2.0) "))
        return [conc_acid, conc_base, volume_acid, list_Vb]
    elif decider == 2:
        conc_base = float(input("\nEnter molarity of Strong base: "))
        conc_acid = float(input("Enter molarity of Strong acid: "))
        volume_base = float(input("\nEnter milliliters of Strong base: "))
        list_Va = str(input("Enter volume of Strong acid ex:(1.0,2.0) "))
        return [conc_base, conc_acid, volume_base, list_Va]

    elif decider == 3:
        conc_acid = float(input("\nEnter molarity of Weak acid : "))
        conc_base = float(input("Enter molarity of Strong base : "))
        ka = float(input("\nEnter ka : "))
        volume_acid = float(input("\nEnter milliliters of Weak acid: "))
        print("Enter volume of Strong base ex:(1.0,2.0) ")
        list_Vb = 10
        return [conc_acid, conc_base, ka, volume_acid, list_Vb]
    elif decider == 4:
        conc_acid = float(input("\nEnter molarity of Weak base : "))
        conc_base = float(input("Enter molarity of Strong acid : "))
        kb = float(input("\nEnter kb : "))
        volume_base = float(input("\nEnter milliliters of Weak base: "))
        list_Va = str(input("Enter volume of Acid base ex:(1.0,2.0) "))
        return [conc_base, conc_acid, kb, volume_base, list_Va]
