from typing import Generator


def median_flow(hydrogram: list) -> tuple:
    hist = {}
    for d in hydrogram:
        bucket = int(d[1])
        hist[bucket] = hist.get(bucket, 0) + 1

    l = sum(hist.values()) / 2
    c = 0
    for b in sorted(hist.keys()):
        if c + hist[b] > l:
            return b, int(hydrogram[c][0])

        c += hist[b]

    return 0, 0

def gen_hydrogram(
    hydrogram: list
) -> Generator[float, None, float]:
    # m, mt = median_flow(hydrogram)
    t = 0
    r = tuple(float(v) for v in hydrogram.pop(0))
    while True:
        # if r[1] <= m and t >= mt:
        #     yield m
        # else:
        if t >= float(hydrogram[0][0]):
            r = tuple(float(v) for v in hydrogram.pop(0))
    
        if len(hydrogram) == 0:
            yield r[1]
            break
    
        yield r[1] + (float(hydrogram[0][1]) - r[1]) / (float(hydrogram[0][0]) - t)
        t += 1

    return 0

