import dataclasses, math

@dataclasses.dataclass(repr=False)
class Error:
    absolute: float
    relative: float

    def __repr__(self) -> str:
        return f"abs: {self.absolute:.2E}, rel: {self.relative:.2E}"


def compute_error(measured: float, expected: float) -> Error:
    absolute = abs(measured - expected)

    if expected != 0:
        relative = absolute / abs(expected)
    elif measured == expected:
        relative = 0
    else:
        relative = float("NaN")

    return Error(absolute, relative)


class AverageError:
    accumulated: Error
    count:       int

    def __init__(self) -> None:
        self.accumulated = Error(0, 0)
        self.count       = 0

    def get(self) -> Error:
        if self.count == 0:
            return Error(0, 0)

        return Error(self.accumulated.absolute / self.count,
                     self.accumulated.relative / self.count)

    def push(self, error: Error) -> None:
        # Do not include nans in the result
        # See: compute_error()
        if math.isnan(error.relative):
            return

        self.accumulated.absolute += error.absolute
        self.accumulated.relative += error.relative

        self.count += 1

    def __repr__(self) -> str:
        return self.get().__repr__()
