import sys

def main():
    if (len(sys.argv) != 2):
        print(len(sys.argv))
        print("This program requires one argument.")
        quit()
    my_name = sys.argv[1]
    print("Hello ",my_name, " !!")
    return


if __name__ == '__main__':

    # run the main function
    main()
