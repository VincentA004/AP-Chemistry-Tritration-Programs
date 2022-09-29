import chemMath as cm 

def main():
    password = "red"
    userPassword = input("Enter password to use Acid-Base Program: ")
    if (password == userPassword):
        check = True
        while (check == True):
            reader = Menu()
            if(reader == 1):
                cm.acid_base_main()
            elif(reader == 2):
                cm.Solution_main()
            elif(reader == 3):
                check = False
        
def Menu():
    print('\n'+"**************************************\n\tACID - BASE PROGRAM\n\tBy: Vincent Allam\n**************************************"+'\n')
    print("1. Basic acid-base equilibrium ")
    print("2. Advanced acid-base equilibrium ")
    print("3. Quit ")
    option_decider = int(input())
    return option_decider

main()