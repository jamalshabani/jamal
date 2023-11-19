import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 simple_obj_static_heat.py -o test1 -tao_monitor -ls 0.055 -lr 5.0 -tao_max_it 5000",
                "python3 simple_obj_static_heat.py -o test2 -tao_monitor -ls 0.06 -lr 5.0 -tao_max_it 5000",
                "python3 simple_obj_static_heat.py -o test3 -tao_monitor -ls 0.065 -lr 5.0 -tao_max_it 5000",
                "python3 simple_obj_static_heat.py -o test4 -tao_monitor -ls 0.07 -lr 5.0 -tao_max_it 5000",
                "python3 simple_obj_static_heat.py -o test5 -tao_monitor -ls 0.075 -lr 5.0 -tao_max_it 5000"]


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
