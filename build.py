import matplotlib.pyplot as plt

data = """
1 - 0.11609359999999999 0.1772274 0.02849292
2 - 0.0631684 0.09418499999999999 0.01689542
4 - 0.0333144 0.0479586 0.00883064
8 - 0.029659 0.0411764 0.0106902
16 - 0.0209674 0.0390352 0.01013908
32 - 0.0277722 0.0580348 0.01552306
64 - 0.078623 0.0802808 0.05710506
128 - 0.1590606 0.1587374 0.12869386
"""

def parse_data(data):
    lines = data.strip().split("\n")
    times = {
        "random-10000": [],
        "corner-10000": [],
        "sparse-50000": []
    }

    for line in lines:
        tokens = line.split()
        times["random-10000"].append(float(tokens[2].strip()))
        times["corner-10000"].append(float(tokens[3].strip()))
        times["sparse-50000"].append(float(tokens[4].strip()))
    
    return times

def calculate_speedup(times):
    speedups = {
        key: [times[key][0] / time for time in value] for key, value in times.items()
    }
    return speedups

def plot_graph(data, title):
    thread_counts = [1, 2, 4, 8, 16, 32, 64, 128]
    for key, values in data.items():
        plt.plot(thread_counts, values, label=key)

    plt.xlabel('Threads')
    plt.ylabel('Speedup')
    plt.title(title)
    plt.xticks(thread_counts)
    plt.legend()
    plt.grid(True)
    plt.savefig(title + ".png")

if __name__ == "__main__":
    times = parse_data(data)
    speedups = calculate_speedup(times)
    plot_graph(speedups, '(simulate) Speedup vs Threads')
