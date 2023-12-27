import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3", "rm -rf test4",
                "python3 circle.py -tao_monitor -tao_max_it 2000 -tao_ls_type armijo -o test1 -k 1.0e-3 -e 1.0e-2 -ls 0.06 -lr 5.0"]


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


# Good
# "python3 circle.py -tao_monitor -tao_max_it 2000 -tao_ls_type armijo -o test1 -k 1.0e-3 -e 1.0e-2 -ls 0.02 -lr 3.0"