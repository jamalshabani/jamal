import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 time_obj_time_heat.py -tao_monitor -tao_max_it 500 -o test1 -k 1.0e-3 -e 1.0e-2 -ls 0.035 -lr 4.0"]


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
