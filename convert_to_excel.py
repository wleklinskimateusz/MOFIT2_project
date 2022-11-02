
import sys


def main():
    # get filename from args
    filename = sys.argv[1]
    with open(filename, 'r') as f:
        with open(f"{filename.replace('.csv', '')}_exel.csv", 'w') as f2:
            f2.write(f.read().replace(',', ';').replace(".", ","))


if __name__ == "__main__":
    main()
