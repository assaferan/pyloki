import numpy as np
from pyloki.utils import transforms

rng = np.random.default_rng(0)

print("Checking: does dropping exactly the TOP order (n_keep = n_params - 1) via")
print("economization differ from naive truncation, per retained-order parity?\n")

for n_params in [3, 4, 5, 6, 7]:
    k_max = n_params - 1
    d_vec = rng.uniform(1, 10, size=n_params)  # [d_kmax, ..., d_1, d_0]
    t_s = 2.3

    naive = d_vec.copy()
    naive[0] = 0.0  # drop the top order directly

    econ = transforms.economize_taylor_params(d_vec, t_s, n_keep=n_params - 1)

    print(f"k_max={k_max} (dropped order {k_max}, parity={'odd' if k_max % 2 else 'even'}):")
    for i in range(n_params):
        order = k_max - i
        diff = econ[i] - naive[i]
        same_parity = (order % 2) == (k_max % 2)
        print(f"  order {order} ({'same' if same_parity else 'diff'} parity as dropped): "
              f"naive={naive[i]:.6f} econ={econ[i]:.6f} diff={diff:.3e}")
    print()
