import os
import glob
import typing

import pandas  as pd
import seaborn as sns


def generate_cpg_style() -> dict:
    BG_COLOR = '#1a1a1a'
    TX_COLOR = '#FFFFFF'

    return {
        'axes.facecolor':    '#121212',
        'axes.edgecolor':    BG_COLOR,
        'axes.labelcolor':   TX_COLOR,
        'text.color':        TX_COLOR,
        'xtick.color':       TX_COLOR,
        'ytick.color':       TX_COLOR,
        'grid.color':        BG_COLOR,
        'figure.facecolor':  BG_COLOR,
        'figure.edgecolor':  BG_COLOR,
        'savefig.facecolor': BG_COLOR,
        'savefig.edgecolor': BG_COLOR,
    }


# pylint: disable=too-many-instance-attributes
class Case:
    def __init__(self, dirpath: str, dt = None, parallel_io: bool = False):
        assert not parallel_io, "Parallel I/O is not supported yet."

        self._dirpath    = dirpath
        self._data       = {}
        self._procs      = set()
        self._timesteps  = set()
        self._timestamps = set()
        self._ndims      = 0
        self._axes       = []
        self._coords     = [set(), set(), set()]
        self._dt         = dt

        self._minmax_time = {}
        self._minmax_step = {}

        for f in glob.glob(os.path.join(self._dirpath, 'D', f'cons.1.*.*.dat')):
            self._procs.add(int(f.split('.')[-3]))
            step = int(f.split('.')[-2])
            self._timesteps.add(step)
            self._timestamps.add(step * (self._dt or 1))

        df_t0_p0 = self._read_csv('cons.1', 0, 0)
        self._ndims = len(df_t0_p0.columns) - 1
        for dim in range(self._ndims):
            self._coords[dim] = set(df_t0_p0.iloc[:, dim])
            self._axes.append(['x', 'y', 'z'][dim])

        for t_step in self._timesteps:
            df = pd.DataFrame()
            for proc in self._procs:
                df = pd.concat([
                    df,
                    self._read_csv(
                        'cons.1', proc, t_step,
                        names=self._axes, usecols=self._axes
                    )
                ])

            self._data[t_step]        = df
            self._minmax_step[t_step] = {}

        for axis in self._axes:
            self._compute_minmax(axis)

    def _read_csv(self, path: str, proc: int, t_step: int, **kwargs) -> pd.DataFrame:
        return pd.read_csv(
            os.path.join(self._dirpath, 'D', f'{path}.{proc:02d}.{t_step:06d}.dat'),
            sep=r'\s+', header=None, **kwargs
        )

    def get_ndims(self)      -> int:  return self._ndims
    def get_coords(self)     -> set:  return self._coords
    def get_timesteps(self)  -> set:  return self._timesteps
    def get_timestamps(self) -> set:  return self._timestamps
    def get_procs(self)      -> set:  return self._procs
    def get_data(self)       -> dict: return self._data

    def define_variable(self, name: str, func: typing.Callable):
        for t_step, data in self._data.items():
            data[name] = data.apply(
                lambda row: func(t_step, row[self._axes]),
                axis=1
            )

        self._compute_minmax(name)

    def load_variable(self, name: str, path: str):
        for t_step in self._timesteps:
            dfs = []
            for proc in self._procs:
                dfs.append(self._read_csv(
                    path, proc, t_step,
                    names=[*self._axes, name]
                ))

            self._data[t_step] = pd.merge(self._data[t_step], pd.concat(dfs))

        self._compute_minmax(name)

    def _compute_minmax(self, varname: str):
        lmins, lmaxs = set(), set()
        for t_step in self._timesteps:
            lmin, lmax = self._data[t_step][varname].min(), self._data[t_step][varname].max()
            self._minmax_step[t_step][varname] = (lmin, lmax)
            lmins.add(lmin); lmaxs.add(lmax)

        self._minmax_time[varname] = (min(lmins), max(lmaxs))

    def get_minmax_time(self, varname: str) -> tuple:
        return self._minmax_time[varname]

    def get_minmax_step(self, varname: str, t_step: int) -> tuple:
        return self._minmax_step[t_step][varname]

    def plot_time(self, varname: str, aggregator: typing.Callable = None, **kwargs):
        if aggregator is None:
            aggregator = lambda x: x.mean()

        return sns.lineplot(
            x=list((self._dt or 1) * t for t in self._timesteps),
            y=[
                aggregator(self._data[t_step][varname])
                for t_step in self._timesteps
            ], **kwargs)

    def plot_step(self, t_step: int, varname: str, axes: str = None, **kwargs):
        axes = axes or self._axes

        if len(axes) == 1:
            return sns.lineplot(self._data[t_step], x=axes[0], y=varname, **kwargs)

        if len(axes) == 2:
            return sns.heatmap(
                self._data[t_step].pivot(
                    index=axes[0], columns=axes[1], values=varname
                ),
                **kwargs
            )

        assert False, "3D plotting is not supported yet."
