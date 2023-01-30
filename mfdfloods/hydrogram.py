from typing import Generator


def gen_hydrogram(
    hydrogram: list[tuple[float, float]]
) -> Generator[float, None, float]:
    t = 0
    r = tuple(float(v) for v in hydrogram.pop(0))
    while True:
        if len(hydrogram) == 0:
            yield r[1]
            break

        if t >= float(hydrogram[0][0]):
            r = tuple(float(v) for v in hydrogram.pop(0))

        yield r[1] + (float(hydrogram[0][1]) - r[1]) / (float(hydrogram[0][0]) - t)
        t += 1

    return 0

