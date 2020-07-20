"""
Hedgehoge Magnetic Skyrmion
Monte Carlo Simulator
"""
import numpy as np
import time
from matplotlib import pyplot as plt
import os


# For periodic boundary condition
class IndexWrapper(object):
    def __init__(self, arr, size):
        if type(size) is not int:
            raise TypeError
        self.size = size
        self.arr = arr
        assert self.arr.ndim == 2

    def modify_tuple(self, t):
        new = ()
        for x in t:
            if type(x) is int:
                new += (x % self.size,)
            elif type(x) is slice:
                assert x.step is None
                x.start %= self.size
                x.stop %= self.size
                new += (x,)
            else:
                new += (x,)
        return new

    def __getitem__(self, item):
        if type(item) is tuple:
            new = self.modify_tuple(item)
            return self.arr.__getitem__(new)
        else:
            return self.arr.__getitem__(item)

    def __setitem__(self, key, value):
        if type(key) is tuple:
            new = self.modify_tuple(key)
            self.arr.__setitem__(new, value)
        else:
            self.arr.__setitem__(key, value)

    def __repr__(self):
        return self.arr.__repr__()


class Parameters(object):
    def __init__(self, j, d, b):
        self.J, self.D, self.B = float(j), float(d), float(b)

    def __repr__(self):
        return "(%.4f, %.4f, %.4f)" % (self.J, self.D, self.B)


class Schedule(object):
    def __init__(self):
        self.dwell = 1000
        self.eps = 1E-1
        self.accepted = 0
        self.tests = 0
        self.T = None
        self.T0 = None

    def init(self, dwell, eps, T0):
        self.dwell = dwell
        self.eps = eps
        self.T = T0
        self.T0 = T0

    def accept_test(self, energy_diff):
        T = self.T
        self.tests += 1
        p = np.exp(-energy_diff * 1.0 / T)
        if energy_diff < 0 or p > np.random.uniform(0.0, 1.0):
            self.accepted += 1
            return True
        else:
            return False

    def __repr__(self):
        return self.__class__.__name__ + str({'dwell': self.dwell, 'eps': self.eps, 'T0': self.T0})

    def update_temp(self):
        pass


class Fast_sa(Schedule):
    def __init__(self):
        super().__init__()
        self.k = 0

    def update_temp(self):
        self.k += 1
        self.T = self.T0 * (0.995 ** self.k)


class DoNothing_sa(Schedule):
    def update_temp(self):
        self.T = self.T0


class Boltzmann_sa(Schedule):
    def __init__(self):
        super().__init__()
        self.k = 0

    def update_temp(self):
        alpha = 0.85
        self.T = self.T0 / (1 + alpha * np.log(1 + self.k))
        self.k += 1


class Hedgehoge:
    class State(object):
        def __init__(self):
            self.cost = None
            self.arr_th = None
            self.arr_ph = None
            self.spin = None
            self.number = None

    def __init__(self, size=24, param=(), options=None):
        if type(size) is not int:
            raise TypeError
        self.size = size
        self.param = Parameters(*param)
        self.state = self.State()

        self.state.arr_th = np.pi * (np.random.rand(self.size, self.size) * 2.0 - 1.0)
        self.state.arr_ph = 2.0 * np.pi * (np.random.rand(self.size, self.size) * 2.0 - 1.0)

        # Schedule
        self.schedule = Boltzmann_sa()
        self.schedule.init(**options)

        self.arr_th = IndexWrapper(self.state.arr_th, self.size)
        self.arr_ph = IndexWrapper(self.state.arr_ph, self.size)

    def __repr__(self):
        return "Hedgehoge Magnetic Skyrmion Simulator\n" \
               + 'Size : ' + str(self.size) + '\n' \
               + repr(self.param)

    @staticmethod
    def o3(th, ph):
        return np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])

    def calc_spin(self):
        th, ph = self.state.arr_th, self.state.arr_ph
        temp = np.ndarray(shape=(self.size, self.size, 3), dtype=np.float32)
        temp[:, :, 0] = np.sin(th) * np.cos(ph)
        temp[:, :, 1] = np.sin(th) * np.sin(ph)
        temp[:, :, 2] = np.cos(th)
        self.state.spin = temp

    def energy(self):
        energy_of_system = 0.0

        arr_th = self.arr_th
        arr_ph = self.arr_ph

        for i in range(self.size):
            for j in range(self.size):
                m00 = self.o3(arr_th[i, j], arr_ph[i, j])
                m10 = self.o3(arr_th[i + 1, j], arr_ph[i + 1, j])
                m01 = self.o3(arr_th[i, j + 1], arr_ph[i, j + 1])

                energy_of_system -= self.param.J * m00 @ m10
                energy_of_system -= self.param.J * m00 @ m01

                energy_of_system += self.param.D * np.cross(m00, m10) @ [0, 1, 0]
                energy_of_system -= self.param.D * np.cross(m00, m01) @ [1, 0, 0]

                energy_of_system -= self.param.B * m00 @ [0, 0, 1]

        return energy_of_system

    def number(self):
        arr_th = self.arr_th
        arr_ph = self.arr_ph

        result = 0.0

        for i in range(self.size):
            for j in range(self.size):
                m00 = self.o3(arr_th[i, j], arr_ph[i, j])
                m10 = self.o3(arr_th[i + 1, j], arr_ph[i + 1, j])
                m01 = self.o3(arr_th[i, j + 1], arr_ph[i, j + 1])
                result += np.dot(m00, np.cross(m10, m01))
        result /= (4.0 * np.pi)
        return result

    def delta(self, i, j, var_th, var_ph):
        arr_th = self.arr_th
        arr_ph = self.arr_ph

        # Original spin
        S_origin = self.o3(arr_th[i, j], arr_ph[i, j])
        # Varied spin
        S_prime = self.o3(arr_th[i, j] + var_th, arr_ph[i, j] + var_ph)

        # Neighborhoods S_r+x, S_r-x, S_r+y, S_r-y
        Sxp = self.o3(arr_th[i + 1, j], arr_ph[i + 1, j])
        Sxn = self.o3(arr_th[i - 1, j], arr_ph[i - 1, j])
        Syp = self.o3(arr_th[i, j + 1], arr_ph[i, j + 1])
        Syn = self.o3(arr_th[i, j - 1], arr_ph[i, j - 1])

        delta_S = S_prime - S_origin

        diff = 0.0
        diff -= self.param.J * np.dot(delta_S, Sxp + Sxn + Syp + Syn)
        diff += self.param.D * np.dot(delta_S, np.cross(Sxp - Sxn, [0, 1, 0]))
        diff -= self.param.D * np.dot(delta_S, np.cross(Syp - Syn, [1, 0, 0]))
        diff -= self.param.B * np.dot(delta_S, [0, 0, 1])

        return diff

    def update(self, i, j, var_th, var_ph):
        self.arr_th[i, j] += var_th
        self.arr_ph[i, j] += var_ph

    def doMC(self):
        for _ in range(self.schedule.dwell):
            for i in range(self.size):
                for j in range(self.size):
                    var_th = self.schedule.eps * np.random.uniform(-1.0, 1.0)
                    var_ph = self.schedule.eps * np.random.uniform(-1.0, 1.0)
                    diff = self.delta(i, j, var_th, var_ph)

                    if self.schedule.accept_test(diff):
                        self.update(i, j, var_th, var_ph)
        self.schedule.update_temp()

        self.state.cost = self.energy()
        self.state.number = self.number()
        # self.calc_spin()

    def save(self, d):
        name = str(self.size) + repr(self.param)
        np.savetxt(d + '/' + 'TH' + name + '.csv', self.state.arr_th, delimiter=',')
        np.savetxt(d + '/' + 'PH' + name + '.csv', self.state.arr_ph, delimiter=',')


if __name__ == '__main__':
    # ======================================================================
    # Lattice size
    s = 24
    # Parameters B D J
    p = (1.0, 6.0 ** 0.5, 2.0)
    # Inner loop, epsilon, Initial temperature
    o = {'dwell': 100, 'eps': 0.2, 'T0': 0.5}
    # ======================================================================
    n = time.strftime("%Y%m%d%H%M%S", time.localtime())

    try:
        os.mkdir(n)
    except FileExistsError:
        pass

    f = open(n + '/' + 'log.txt', 'wt')
    lattice = Hedgehoge(size=s,
                        param=p,
                        options=o)
    print(lattice)
    print(lattice, file=f)
    print(lattice.schedule)
    print(lattice.schedule, file=f)

    f.close()

    try:
        # Outer loop
        max_iter = 100

        my_dpi = 96
        plt.ion()
        fig = plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
        plt.pause(5.0)

        # Monte Carlo simulation
        energy_log = []
        for it in range(max_iter):
            t_n = time.time()

            lattice.doMC()

            t_n = time.time() - t_n
            print('Iter :', it, 'Energy : %.4f' % lattice.state.cost,
                  'Number : %.4f' % lattice.state.number,
                  'Temp : %.4f' % lattice.schedule.T,
                  'Time : %.4f' % t_n)

            lattice.save(n)

            # spin = lattice.state.spin

            energy_log.append(lattice.state.cost)

            arr_th_n = lattice.arr_th
            arr_ph_n = lattice.arr_ph
            d_arr = np.zeros((s, s))
            for ix in range(d_arr.shape[0]):
                for jx in range(d_arr.shape[1]):
                    m00_n = Hedgehoge.o3(arr_th_n[ix, jx], arr_ph_n[ix, jx])
                    m10_n = Hedgehoge.o3(arr_th_n[ix + 1, jx], arr_ph_n[ix + 1, jx])
                    m01_n = Hedgehoge.o3(arr_th_n[ix, jx + 1], arr_ph_n[ix, jx + 1])
                    d_arr[ix, jx] = np.dot(m00_n, np.cross(m01_n, m10_n))
            plt.gca().clear()
            ax1.pcolormesh(d_arr.T, cmap='RdBu', vmin=-1.0, vmax=1.0, shading='gouraud')  # shading='gouraud'
            ax1.set_aspect('equal')
            ax2.plot(energy_log)
            plt.pause(1E-5)
            fig.savefig(n + '/' + 'Result.png', dpi=my_dpi)
    except KeyboardInterrupt:
        pass
