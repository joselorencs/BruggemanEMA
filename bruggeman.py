"""Calculate effective n and k of a mixed layer using Bruggeman.

The program opens a dialog box to select one or two TXT files. Each TXT file must have a header and three numerical columns: wavelength (nm), n, and k.

If only one file is selected, the second material is assumed to be air: n = 1 and k = 0 across the entire spectral range.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path
from tkinter import Tk, messagebox, simpledialog
from tkinter.filedialog import askopenfilenames, asksaveasfilename

import numpy as np


DEPOLARIZATION_Q = 0.5
DEFAULT_F2 = 0.5


def read_material_txt(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read wavelength, n, and k from a TXT file with one header row."""
    try:
        text = path.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError:
        text = path.read_text(encoding="latin-1")

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"{path.name} does not contain enough header and data rows.")

    data_lines = lines[1:]
    delimiter = detect_delimiter(data_lines[0])
    data = np.genfromtxt(io.StringIO("\n".join(data_lines)), delimiter=delimiter)

    if data.ndim == 1:
        if data.size != 3:
            raise ValueError(f"{path.name} must have exactly 3 columns.")
        data = data.reshape(1, 3)

    if data.shape[1] < 3:
        raise ValueError(f"{path.name} must have 3 columns: wavelength (nm), n, and k.")

    data = data[:, :3]
    if not np.all(np.isfinite(data)):
        raise ValueError(f"{path.name} contains non-numeric or infinite values.")

    order = np.argsort(data[:, 0])
    data = data[order]

    wavelength = data[:, 0]
    if np.any(np.diff(wavelength) == 0):
        raise ValueError(f"{path.name} contains duplicated wavelength values.")

    return wavelength, data[:, 1], data[:, 2]


def detect_delimiter(first_data_line: str) -> str | None:
    if ";" in first_data_line:
        return ";"
    if "," in first_data_line:
        return ","
    return None


def align_materials(
    material_1: tuple[np.ndarray, np.ndarray, np.ndarray],
    material_2: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, bool]:
    w1, n1, k1 = material_1

    if material_2 is None:
        return w1, n1, k1, np.ones_like(w1), np.zeros_like(w1), False

    w2, n2, k2 = material_2
    same_grid = len(w1) == len(w2) and np.allclose(w1, w2, rtol=1e-9, atol=1e-9)
    if same_grid:
        return w1, n1, k1, n2, k2, False

    overlap_min = max(float(w1.min()), float(w2.min()))
    overlap_max = min(float(w1.max()), float(w2.max()))
    if overlap_min >= overlap_max:
        raise ValueError("The two materials do not have a common spectral range.")

    mask = (w1 >= overlap_min) & (w1 <= overlap_max)
    if not np.any(mask):
        raise ValueError("There are no material 1 points within the common range.")

    wavelength = w1[mask]
    n2_interp = np.interp(wavelength, w2, n2)
    k2_interp = np.interp(wavelength, w2, k2)

    return wavelength, n1[mask], k1[mask], n2_interp, k2_interp, True


def bruggeman_effective_epsilon(
    epsilon_1: np.ndarray,
    epsilon_2: np.ndarray,
    f2: float,
    q: float = DEPOLARIZATION_Q,
) -> np.ndarray:
    """Solve the symmetric Bruggeman equation for the effective epsilon."""
    f1 = 1.0 - f2
    a = q
    b = 1.0 - q

    if np.isclose(b, 0.0):
        raise ValueError("q cannot be 1 for this equation form.")

    beta = f1 * (b * epsilon_1 - a * epsilon_2) + f2 * (b * epsilon_2 - a * epsilon_1)
    discriminant = beta * beta + 4.0 * a * b * epsilon_1 * epsilon_2
    root = np.sqrt(discriminant)

    solution_plus = (beta + root) / (2.0 * b)
    solution_minus = (beta - root) / (2.0 * b)

    # The physical branch is usually close to the linear permittivity average.
    reference = f1 * epsilon_1 + f2 * epsilon_2
    use_plus = np.abs(solution_plus - reference) <= np.abs(solution_minus - reference)
    return np.where(use_plus, solution_plus, solution_minus)


def nk_to_epsilon(n: np.ndarray, k: np.ndarray) -> np.ndarray:
    return (n + 1j * k) ** 2


def epsilon_to_nk(epsilon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    epsilon_array = np.asarray(epsilon, dtype=complex)
    flat_epsilon = epsilon_array.ravel()
    refractive_index = np.empty_like(flat_epsilon, dtype=complex)
    previous_value: complex | None = None
    tolerance = 1e-10

    for index, epsilon_value in enumerate(flat_epsilon):
        root = np.sqrt(epsilon_value)
        candidates = np.array([root, -root], dtype=complex)
        k_values = np.imag(candidates)
        passive_candidates = np.where(k_values >= -tolerance)[0]

        if len(passive_candidates) == 1 and abs(k_values[passive_candidates[0]]) > tolerance:
            selected = candidates[passive_candidates[0]]
        else:
            if previous_value is None:
                nonnegative_n = np.where(np.real(candidates) >= -tolerance)[0]
                if len(nonnegative_n) > 0:
                    selected = candidates[nonnegative_n[0]]
                elif len(passive_candidates) > 0:
                    selected = candidates[passive_candidates[0]]
                else:
                    selected = candidates[np.argmax(k_values)]
            else:
                if len(passive_candidates) > 0:
                    candidate_pool = candidates[passive_candidates]
                else:
                    candidate_pool = candidates
                selected = candidate_pool[np.argmin(np.abs(candidate_pool - previous_value))]

        if np.imag(selected) < 0:
            if np.imag(selected) >= -tolerance:
                selected = complex(np.real(selected), 0.0)
            else:
                selected = -selected

        refractive_index[index] = selected
        previous_value = selected

    refractive_index = refractive_index.reshape(epsilon_array.shape)
    n = np.real(refractive_index)
    k = np.maximum(np.imag(refractive_index), 0.0)

    n = np.where(np.abs(n) < 1e-12, 0.0, n)
    k = np.where(np.abs(k) < 1e-12, 0.0, k)
    return n, k


def save_output(path: Path, wavelength: np.ndarray, n: np.ndarray, k: np.ndarray) -> None:
    output = np.column_stack((wavelength, n, k))
    np.savetxt(
        path,
        output,
        fmt="%.10g",
        delimiter="\t",
        header="wavelength (nm)\tn\tk",
        comments="",
    )


def run() -> None:
    root = Tk()
    root.withdraw()
    root.update()

    try:
        selected = askopenfilenames(
            title="Select 1 or 2 material TXT files",
            filetypes=[("TXT files", "*.txt"), ("All files", "*.*")],
        )

        paths = [Path(path) for path in selected]
        if not paths:
            return
        if len(paths) > 2:
            messagebox.showerror("Invalid selection", "Select only 1 or 2 TXT files.")
            return

        material_1 = read_material_txt(paths[0])
        material_2 = read_material_txt(paths[1]) if len(paths) == 2 else None

        f2 = simpledialog.askfloat(
            "Volume fraction",
            "Volume fraction of material 2\n(air if you selected only one TXT file):",
            initialvalue=DEFAULT_F2,
            minvalue=0.0,
            maxvalue=1.0,
            parent=root,
        )
        if f2 is None:
            return

        q = simpledialog.askfloat(
            "Depolarization factor",
            "Depolarization factor q:",
            initialvalue=DEPOLARIZATION_Q,
            minvalue=0.0,
            maxvalue=1.0,
            parent=root,
        )
        if q is None:
            return

        wavelength, n1, k1, n2, k2, interpolated = align_materials(material_1, material_2)

        epsilon_1 = nk_to_epsilon(n1, k1)
        epsilon_2 = nk_to_epsilon(n2, k2)
        epsilon_eff = bruggeman_effective_epsilon(epsilon_1, epsilon_2, f2, q)
        n_eff, k_eff = epsilon_to_nk(epsilon_eff)

        default_name = f"{paths[0].stem}_bruggeman_mix.txt"
        output_path = asksaveasfilename(
            title="Save mixed-layer TXT file",
            initialdir=str(paths[0].parent),
            initialfile=default_name,
            defaultextension=".txt",
            filetypes=[("TXT files", "*.txt"), ("All files", "*.*")],
        )
        if not output_path:
            return

        save_output(Path(output_path), wavelength, n_eff, k_eff)

        note = ""
        if len(paths) == 1:
            note = "\nMaterial 2: air (n=1, k=0)."
        elif interpolated:
            note = "\nMaterial 2 was interpolated to the common range of material 1."

        messagebox.showinfo(
            "Process completed",
            f"File saved successfully:\n{output_path}{note}",
        )
    except Exception as exc:
        messagebox.showerror("Error", str(exc))
        raise
    finally:
        root.destroy()


if __name__ == "__main__":
    try:
        run()
    except Exception:
        sys.exit(1)
