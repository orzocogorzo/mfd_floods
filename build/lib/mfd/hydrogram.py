def hydrogram (break_flow, base_flow, break_seconds):
    # for t in range(break_seconds):
    #     yield break_flow * break_flow**(-4*t/(break_seconds*0.3)) + base_flow
    t = 0
    while True:
        yield break_flow * break_flow ** (-1 * t / (break_seconds * 0.3)) + base_flow
        t += 1