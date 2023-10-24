import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 simple_obj_time_heat.py -tao_type bncg -tao_max_funcs 10000 -tao_monitor -tao_max_it 1000 -tao_ls_type more-thuente -m 'trajectory.msh' -o 'test1' -er 1.0e1 -es 1.0e-1 -lr 5.0 -ls 0.05 -vr 0.3 -vs 0.3 -k 2.5e-3 -e 3.0e-2 -p 2.0 -q 1.0",
                "python3 simple_obj_time_heat.py -tao_type bncg -tao_max_funcs 10000 -tao_monitor -tao_max_it 1000 -tao_ls_type more-thuente -m 'trajectory.msh' -o 'test2' -er 1.0e1 -es 1.0e-1 -lr 4.0 -ls 0.05 -vr 0.3 -vs 0.3 -k 2.5e-3 -e 3.0e-2 -p 2.0 -q 1.0"]


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
