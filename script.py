import subprocess
import re
from statistics import mean
import matplotlib.pyplot as plt


def run_script(threads,cmd):
    cmd = f"OMP_NUM_THREADS={threads} " + cmd
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)  # 5 minutes timeout
    except subprocess.TimeoutExpired:
        print(f"Command timed out after 5 minutes when running with {threads} threads!")
        return None
    return result.stdout

thread_counts = [1, 2, 4, 8, 16, 32, 64, 128]

def parse_build_speedup(output):

    lines = output.split("\n")
    tree = []
    simulate = []
    for line in lines:
        if "iteration" in line:
            tree.append(float(line.split()[-3][:-2]))
            simulate.append(float(line.split()[-1][:-1]))
            
    tree = mean(tree)
    simulate = mean(simulate)
    return tree, simulate
            
def plot_graph(data, title):
    for key, values in data.items():
        plt.plot(thread_counts, values, label=key)
    
    plt.figure()
    plt.xlabel('Threads')
    plt.ylabel('Speedup')
    plt.title(title)
    plt.xticks(thread_counts)
    plt.legend()
    plt.grid(True)
    plt.savefig(title + ".png")


def main():
    cmds = [f"./nbody-release -task -n 10000 -i 5 -in ./src/benchmark-files/random-10000-init.txt -s 100 -parallelfor"
,f"./nbody-release -task -n 10000 -i 5 -in ./src/benchmark-files/corner-10000-init.txt -s 100 -parallelfor",
f"./nbody-release -task -n 10000 -i 50 -in ./src/benchmark-files/sparse-50000-init.txt -s 5 -parallelfor"]

    print("------------")
    print("Threads\tAverage Build Speedup")
    for threads in thread_counts:
        build = []
        sims = []
        
        for cmd in cmds:
            output = run_script(threads,cmd)
            tree, sim = parse_build_speedup(output)
            build.append(tree)
            sims.append(sim)
        # add to csv file
        
        with open('build_speedup.csv', 'a') as f:
            f.write(f"{threads} - {' '.join(map(str,build))}\n") 
        with open('total_speedup.csv', 'a') as f:
            f.write(f"{threads} - {' '.join(map(str,sims))}\n") 

if __name__ == "__main__":
    main()
