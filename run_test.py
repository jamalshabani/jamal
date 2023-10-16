
import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3", "rm -rf test4", "rm -rf test5",
                "python3 inverse_problem.py -n 500 -o 'test1' -er 1.0 -es 1.0e-2 -vr 0.3 -vs 0.3 -k 1.0e-4 -e 4.0e-3 -p 2.0 -q 1.0"]



i = 1
for program in program_list:
    print("------------------------------------------------------------------------------")
    print("")
    print("Running test #{}".format(i))
    print("")
    print(program)
    print("")
    subprocess.run(program, shell = True)
    i = i + 1
