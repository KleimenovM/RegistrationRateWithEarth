import numpy as np
import matplotlib.pyplot as plt

from transmission_function import TransmissionFunction


def check_transmission():

    tf = TransmissionFunction()

    angle_numbers = [1, 60, 120, 179]

    step = 15

    tf_matrix_emu = np.log10(tf.table_data[1] + 1e-15)[:, ::step, ::step]
    tf_matrix_tau = np.log10(tf.table_data[2] + 1e-15)[:, ::step, ::step]

    size = tf_matrix_emu.shape[-1]

    e = []
    cur_value, new_value = 0, 0
    for i, x in enumerate(tf.lg_energy[::step]):
        new_value = int(x)
        if new_value != cur_value:
            e.append(new_value)
        else:
            e.append('')
        cur_value = new_value

    plt.figure(figsize=(8, 7))

    for i in range(len(angle_numbers)):
        plt.subplot(2, 2, i+1)
        plt.imshow(tf_matrix_emu[angle_numbers[i]].T, vmin=-15, vmax=2, origin='lower')
        plt.xticks(np.arange(size), e)
        plt.yticks(np.arange(size), e)
        if i in [0, 2]:
            plt.ylabel(r"$\lg (E_{out} / GeV)$", fontsize=12)
        if i in [2, 3]:
            plt.xlabel(r"$\lg (E_{in} / GeV)$", fontsize=12)

        plt.tick_params(labelsize=12)

        if i % 2 == 1:
            plt.colorbar()

    plt.tight_layout()
    plt.show()

    return


if __name__ == '__main__':
    check_transmission()

