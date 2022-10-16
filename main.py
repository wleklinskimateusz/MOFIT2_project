from utils import global_index_node, print_node


L = 100  # nm
N = 2


def main():
    for element_idx in range(4*N**2):
        print(f"[ELEMENT]: {element_idx}")
        for local_node in range(4):
            global_node = global_index_node(local_node, N, element_idx)
            print_node(global_node, N, L, element_idx)
            print()


if __name__ == "__main__":
    main()
