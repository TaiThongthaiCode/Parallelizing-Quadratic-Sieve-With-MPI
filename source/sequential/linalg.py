"""
"Candice" - Abraham Lincoln
"""

import qs_helper as qs

def main():
    N = 121651
    bit_matrix = []
    smooth_nums = []

    with open("Expo_Matrix.txt") as f:
        for line in f:
            relation = [int(line[i]) for i in range(len(line)-1)]
            bit_matrix.append(relation)
        f.close()

    with open("Smooth_Num.txt") as f:
        for line in f:
            smooth = int(line)
            smooth_nums.append(smooth)

        f.close()

    # with open("factorbase.txt") as f:
    #     for line in f:
    #         factor_base.append(int(line.split()[0]))

    #     f.close()

    # bit_matrix = qs.build_matrix(smooth_nums, factor_base)
    print("SIZE BIT", len(bit_matrix))
    print("SIZE SMOOTH", len(smooth_nums))



    qs.transpose(bit_matrix)
    bit_matrix.reverse()
    smooth_nums.reverse()
    # print(smooth_nums)


    a, b = qs.lin_alg_solve(bit_matrix, smooth_nums, N)

    print(a)
    print(b)

    return


if __name__ == "__main__":
    main()
