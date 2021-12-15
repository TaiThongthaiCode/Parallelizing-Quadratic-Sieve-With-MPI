def main():

    smooth_nums = []
    factorbase = []
    expo_matrix = []


    with open("Smooth_Num.txt") as sm_file:
        while True:
            line = sm_file.readline()
            if line == "":
                break
            smooth_nums.append(int(line))
    sm_file.close()

    with open("factorbase.txt") as fb_file:
        while True:
            line = fb_file.readline()
            if line == "":
                break
            factorbase.append(int(line.split()[0]))
    fb_file.close


    with open("Power_Matrix.txt") as em_file:
        while True:
            temp = []
            line = em_file.readline()
            if line == "":
                break
            for power in line:
                if power != "\n":
                    temp.append(int(power))
            expo_matrix.append(temp)


    for i in range(len(smooth_nums)):
        res = check(smooth_nums[i], expo_matrix[i], factorbase)
        if not res:
            print("Error at", smooth_nums[i])
    print("All good!")

def check(smooth, relation, factorbase):

    result = 1
    for i in range(len(relation)):
        result *= factorbase[i]**relation[i]

    if result == smooth:
        return True
    else:
        return False






if __name__ == "__main__":
    main()
