import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import tifffile
from sklearn.linear_model import LassoLars


# -------------------------------------------------------------------------
#  Nonnegative LASSO (ADMM) – MATLAB-style using LARS
# -------------------------------------------------------------------------
def nneg_lasso(y, S, lam, a=1.0, n_iter=5):
    """
    Lasso spectral unmixing with non-negativity constraint (ADMM),
    Python version of your MATLAB nneg_lasso.

    Parameters
    ----------
    y : ndarray, shape (Nx, Ny, Nz)
        Hyperspectral image stack.
    S : ndarray, shape (Nz, k)
        Pure component spectra (columns = components).
    lam : float
        Sparsity parameter (same role as MATLAB 'lambda').
    a : float
        ADMM parameter controlling convergence speed.
    n_iter : int
        Number of ADMM iterations.

    Returns
    -------
    C : ndarray, shape (Nx, Ny, k)
        Concentration maps.
    """
    y = np.asarray(y, dtype=float)
    Nx, Ny, Nz = y.shape
    _, k = S.shape

    # Main variables
    C = np.zeros((Nx, Ny, k), dtype=float)
    u = np.zeros_like(C)
    vhat = np.zeros_like(C)
    R_positive = np.zeros_like(C)

    N = Nx * Ny * Nz

    # Build S_tilde as in MATLAB: [S; eye(k)]
    S_tilde = np.vstack([S, np.eye(k)])        # (Nz + k) x k
    n_samples = S_tilde.shape[0]               # Nz + k

    # Standardize predictors like MATLAB (Standardize = true)
    X = S_tilde
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0, ddof=0)
    X_std[X_std == 0] = 1.0                   # avoid divide by zero
    X_scaled = (X - X_mean) / X_std

    # LARS alpha: objective 0.5||y - Xb||^2 + alpha ||b||_1
    # MATLAB: (1/(2n))||y - Xb||^2 + lambda ||b||_1
    # equivalent if alpha = n * lambda
    alpha = lam * n_samples

    print("Iter \t residualC \t residualv \t residualu")
    for it in range(n_iter):
        C_old = C.copy()
        vhat_old = vhat.copy()
        u_old = u.copy()

        ctilde = vhat - u

        # --- Update C (pixel-wise LASSO) ---
        for i in range(Nx):
            for j in range(Ny):
                # Same logic as MATLAB: solve only when needed
                if it == 0 or np.min(C_old[i, j, :]) < 0:
                    y_sp = y[i, j, :]  # (Nz,)
                    rhs = np.concatenate(
                        [y_sp, np.sqrt(a) * ctilde[i, j, :]]
                    )                   # (Nz + k,)

                    # Center response (MATLAB lasso subtracts mean of y)
                    rhs_centered = rhs - rhs.mean()

                    # LARS-based LASSO, no extra normalization or intercept
                    model = LassoLars(
                        alpha=alpha,
                        fit_intercept=False,
                        normalize=False
                    )
                    model.fit(X_scaled, rhs_centered)

                    # Coeffs are for standardized X; undo scaling
                    beta_std = model.coef_          # shape (k,)
                    beta = beta_std / X_std         # back to original scale

                    C[i, j, :] = beta

        # --- Update vhat (non-negativity projection) ---
        vhat = np.maximum(C + u, R_positive)

        # --- Update u ---
        u = u + (C - vhat)

        # Residuals (same formulas as MATLAB)
        residualC = np.linalg.norm(C - C_old) / np.sqrt(N)
        residualv = np.linalg.norm(vhat - vhat_old) / np.sqrt(N)
        residualu = np.linalg.norm(u - u_old) / np.sqrt(N)

        print(
            f"{it+1:3d} \t {residualC:3.5e} \t {residualv:3.5e} \t {residualu:3.5e}"
        )

    return C


# -------------------------------------------------------------------------
#  Main script
# -------------------------------------------------------------------------
def main():
    # ---------------- Load pure chemicals ----------------
    rpf = (2994 - 2913) / (76 - 40)  # Use DMSO to calibrate
    Raman_shift = np.linspace(2913 - 40 * rpf, 2994 + 22 * rpf, 100)


    # Adjust keys if different in your .txt file
    BSA = np.loadtxt('BSA_ref.txt').squeeze()
    TAG = np.loadtxt('TAG_ref.txt').squeeze()
    CHL = np.loadtxt('CHL_ref.txt').squeeze()
    RNA = np.loadtxt('RNA_ref.txt').squeeze()
    GLU = np.loadtxt('GLU_ref.txt').squeeze()

    def normalize(x):
        x = x.astype(float)
        return (x - x.min()) / (x.max() - x.min())

    BSA_n = normalize(BSA)
    TAG_n = normalize(TAG)
    CHL_n = normalize(CHL)
    RNA_n = normalize(RNA)
    GLU_n = normalize(GLU)

    # Quick check plot (optional)
    plt.figure()
    plt.plot(Raman_shift, BSA_n, linewidth=1)
    plt.plot(Raman_shift, TAG_n + 1, linewidth=1)
    plt.plot(Raman_shift, CHL_n + 2, linewidth=1)
    plt.plot(Raman_shift, RNA_n + 3, linewidth=1)
    plt.plot(Raman_shift, GLU_n + 4, linewidth=1)
    plt.xlim([2820, 3050])
    ax = plt.gca()
    ax.minorticks_on()
    plt.legend(["BSA", "TAG", "CHL", "RNA", "GLU"])
    plt.xlabel(r"Raman shift (cm$^{-1}$)")
    plt.ylabel("Int (a.u.)")
    plt.tight_layout()
    plt.show()

    # ---------------- Load target TIFF stack ----------------
    directory = "./"
    target_img = "U87_CHL_3_CH_raw.tif"
    wavelength = 100  # Nz

    stack = tifffile.imread(os.path.join(directory, target_img))

    # tifffile usually gives (Nz, Nx, Ny) for multipage
    if stack.ndim != 3:
        raise ValueError("Expected a 3D TIFF stack (multipage).")

    if stack.shape[0] == wavelength:
        # (Nz, Nx, Ny) -> (Nx, Ny, Nz) to match MATLAB
        y = np.transpose(stack, (1, 2, 0))
    elif stack.shape[2] == wavelength:
        y = stack
    else:
        raise ValueError("Unexpected TIFF stack shape for given wavelength.")

    Nz = wavelength
    Nx, Ny, _ = y.shape

    # ---------------- Background subtraction ----------------
    y_sum = np.mean(y, axis=2)
    plt.figure()
    plt.hist(y_sum.ravel(), bins=50)
    plt.title("Histogram of mean intensity")
    plt.show()

    BGmask = np.zeros_like(y_sum, dtype=bool)
    BG = y_sum < 0.30
    BGmask[BG] = True

    plt.figure()
    plt.imshow(BGmask, cmap="gray")
    plt.title("Background Mask")
    plt.axis("off")
    plt.show()

    BG_spectrum = np.zeros(Nz, dtype=float)
    for i in range(Nz):
        y_temp = y[:, :, i]
        BG_spectrum[i] = y_temp[BGmask].mean()

    plt.figure()
    plt.plot(BG_spectrum)
    plt.title("Background Spectrum")
    plt.xlabel("Channel")
    plt.ylabel("Mean background")
    plt.show()

    # (normalized background spectrum is kept just in case, not needed below)
    _ = normalize(BG_spectrum)

    # subtract scalar BG_spectrum(i) slice-by-slice, like MATLAB
    y_sub = np.empty_like(y, dtype=float)
    for i in range(Nz):
        y_sub[:, :, i] = y[:, :, i] - BG_spectrum[i]

    # Shift to be non-negative
    y_sub = y_sub - y_sub.min()

    # ---------------- Lambda values ----------------
    l_BSA = 1e-4
    l_TAG = 1e-3
    l_CHL = 2.5e-3
    l_RNA = 2.5e-1
    l_GLU = 15e-2

    lambdas = [l_BSA, l_TAG, l_CHL, l_RNA, l_GLU]
    update = [1, 1, 1, 1, 1]  # which component to solve

    # Reference spectra S: (Nz, k)
    ref = np.column_stack([BSA_n, TAG_n, CHL_n, RNA_n, GLU_n])

    a = 1.0        # ADMM parameter
    iter_admm = 5  # # of ADMM iterations

    k = 3
    C = np.zeros((Nx, Ny, k), dtype=float)

    # ---------------- Nonnegative LASSO unmixing ----------------
    print("\n=== Starting nonnegative LASSO unmixing ===")
    # NOTE: nneg_lasso internally solves for all k components each call,
    # but we call it once per lambda to mirror your MATLAB structure.
    for ll in range(k):
        if update[ll] == 1:
            L = lambdas[ll]
            print(f"\nRunning nneg_lasso for component {ll+1} with lambda={L}")
            C_full = nneg_lasso(y_sub, ref, L, a=a, n_iter=iter_admm)
            C[:, :, ll] = C_full[:, :, ll]

    # ---------------- Display results ----------------
    c1 = C[:, :, 0].ravel()
    disp_min = np.percentile(c1, 0.3)
    disp_max = np.percentile(c1, 99.7)
    clims = (disp_min, disp_max)

    titles = ["BSA", "TAG", "CHL", "RNA", "GLU"]

    # Global check
    print("C global min:", np.nanmin(C), "max:", np.nanmax(C))

    for idx in range(C.shape[2]):
        img = C[:, :, idx]

        # Skip if all zeros / NaNs
        print(
            f"{titles[idx]} map: min={np.nanmin(img)}, max={np.nanmax(img)}"
        )

        plt.figure()
        # Use per-map percentiles for contrast
        vmin = np.percentile(img, 1)
        vmax = np.percentile(img, 99)
        if vmin == vmax:
            # fallback: let imshow autoscale
            plt.imshow(img, cmap="gray")
        else:
            plt.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)

        plt.title(titles[idx])
        plt.axis("off")
        plt.colorbar(label="concentration (a.u.)")
        plt.tight_layout()
        plt.show()

    # ---------------- Save TXT chemical maps ----------------
    opt_filepath = "chemical_maps_CH_v2"
    os.makedirs(opt_filepath, exist_ok=True)

    filename = target_img
    ext = ".txt"

    BSA_map = C[:, :, 0]
    TAG_map = C[:, :, 1]
    CHL_map = C[:, :, 2]
    # RNA_map = C[:, :, 3]
    # GLU_map = C[:, :, 4]

    np.savetxt(
        os.path.join(opt_filepath, f"{filename}_BSA{ext}"),
        BSA_map,
        delimiter="\t"
    )
    np.savetxt(
        os.path.join(opt_filepath, f"{filename}_TAG{ext}"),
        TAG_map,
        delimiter="\t"
    )
    np.savetxt(
        os.path.join(opt_filepath, f"{filename}_CHL{ext}"),
        CHL_map,
        delimiter="\t"
    )
    # np.savetxt(
    #     os.path.join(opt_filepath, f"{filename}_RNA{ext}"),
    #     RNA_map,
    #     delimiter="\t"
    # )
    # np.savetxt(
    #     os.path.join(opt_filepath, f"{filename}_GLU{ext}"),
    #     GLU_map,
    #     delimiter="\t"
    # )

    print("Saved chemical maps to:", opt_filepath)


if __name__ == "__main__":
    main()
