# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "polars",
#     "matplotlib",
#     "scienceplots",
#     "numpy",
# ]
# ///
import polars as pl
import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Set plotting style
try:
    import scienceplots
    plt.style.use(['science', 'ieee'])
    plt.rcParams.update({
        'text.usetex': False,
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif', 'serif'],
    })
except ImportError:
    plt.style.use('classic')

def plot_m_vs_csi():
    results_dir = "results"
    csv_file = os.path.join(results_dir, "results_aggregate.csv")
    if not os.path.exists(csv_file): return

    try:
        q = pl.scan_csv(csv_file)
        df = q.with_columns([
            pl.col("b").cast(pl.Utf8).str.replace("d", "e", literal=True).cast(pl.Float64, strict=False).alias("b_val"),
            pl.col("csi").cast(pl.Utf8).str.replace("d", "e", literal=True).cast(pl.Float64, strict=False).alias("csi_val"),
            pl.col("mmax").cast(pl.Float64),
            pl.col("rmax").cast(pl.Float64)
        ]).collect()
    except Exception: return

    if df.height == 0: return

    models = df["model"].unique().to_list()
    
    for model in models:
        mod_df = df.filter(pl.col("model") == model)
        b_vals = mod_df["b_val"].unique().sort().to_list()
        
        fig, ax = plt.subplots(figsize=(6, 4))
        
        for j, b in enumerate(b_vals):
            subset = mod_df.filter(pl.col("b_val") == b).sort("csi_val")
            x = subset["csi_val"].to_list()
            y = subset["mmax"].to_list()
            
            exp = int(np.log10(b)) if b > 0 else 0
            label_str = fr"$B=10^{{{exp}}}$ G"
            ax.plot(x, y, 
                    linestyle='None', marker='o', markersize=3,
                    label=label_str)

        ax.set_title(f"Model: {model} - Max Mass")
        ax.set_xscale('log')
        ax.set_xlabel(r"log$_{10}(\xi)$ [fm]")
        ax.set_ylabel(r"$M_{max}$ [$M_\odot$]")
        ax.legend(fontsize='small', title="Magnetic Field")
        ax.grid(True, alpha=0.3)
        
        out_root = f"mmax_evolution_{model}"
        pdf_path = os.path.join(results_dir, f'{out_root}.pdf')
        png_path = os.path.join(results_dir, f'{out_root}.png')
        
        plt.tight_layout()
        plt.savefig(pdf_path)
        # plt.savefig(png_path, dpi=300)
        print(f"Saved: {pdf_path}")
        plt.close(fig)

def plot_m_vs_rmax():
    results_dir = "results"
    csv_file = os.path.join(results_dir, "results_aggregate.csv")
    if not os.path.exists(csv_file): return

    try:
        q = pl.scan_csv(csv_file)
        df = q.with_columns([
            pl.col("b").cast(pl.Utf8).str.replace("d", "e", literal=True).cast(pl.Float64, strict=False).alias("b_val"),
            pl.col("csi").cast(pl.Utf8).str.replace("d", "e", literal=True).cast(pl.Float64, strict=False).alias("csi_val"),
            pl.col("mmax").cast(pl.Float64),
            pl.col("rmax").cast(pl.Float64)
        ]).collect()
    except Exception: return

    if df.height == 0: return

    models = df["model"].unique().to_list()
    
    for model in models:
        mod_df = df.filter(pl.col("model") == model)
        b_vals = mod_df["b_val"].unique().sort().to_list()
        
        fig, ax = plt.subplots(figsize=(6, 4))
        
        for j, b in enumerate(b_vals):
            subset = mod_df.filter(pl.col("b_val") == b).sort("csi_val")
            x = subset["csi_val"].to_list()
            y = subset["rmax"].to_list()
            
            exp = int(np.log10(b)) if b > 0 else 0
            label_str = fr"$B=10^{{{exp}}}$ G"
            ax.plot(x, y, 
                    linestyle='None', marker='o', markersize=3,
                    label=label_str)

        ax.set_title(f"Model: {model} - Radius at Max Mass")
        ax.set_xscale('log')
        ax.set_xlabel(r"log$_{10}(\xi)$ [fm]")
        ax.set_ylabel(r"$R(M_{max})$ [km]")
        ax.legend(fontsize='small', title="Magnetic Field")
        ax.grid(True, alpha=0.3)
        
        out_root = f"rmax_evolution_{model}"
        pdf_path = os.path.join(results_dir, f'{out_root}.pdf')
        png_path = os.path.join(results_dir, f'{out_root}.png')
        
        plt.tight_layout()
        plt.savefig(pdf_path)
        # plt.savefig(png_path, dpi=300)
        print(f"Saved: {pdf_path}")
        plt.close(fig)

def plot_mr_curves_with_inset():
    results_dir = "results"
    output_dir = "output"
    csv_file = os.path.join(results_dir, "results_aggregate.csv")
    if not os.path.exists(csv_file): return

    try:
        df = pl.read_csv(csv_file)
    except Exception: return

    if df.height == 0: return

    models = df["model"].unique().to_list()

    for model in models:
        mod_df = df.filter(pl.col("model") == model)
        b_list = mod_df["b"].unique().to_list()
        
        for b_str in b_list:
            # Prepare Label
            b_clean = b_str.replace("d", "e")
            try:
                b_val = float(b_clean)
                exp = int(np.log10(b_val)) if b_val > 0 else 0
                b_latex = fr"10^{{{exp}}}"
            except:
                b_latex = b_str

            # Filter csi list (strings)
            csi_list = mod_df.filter(pl.col("b") == b_str)["csi"].to_list()
            
            # Downsample if too many
            if len(csi_list) > 5:
                indices = np.linspace(0, len(csi_list)-1, 5, dtype=int)
                selected_csi = [csi_list[i] for i in indices]
            else:
                selected_csi = csi_list

            fig, ax = plt.subplots(figsize=(3.5, 4))
            inset_curves = []
            
            # Setup colors
            cm = plt.get_cmap('viridis')
            colors = [cm(1.*i/len(selected_csi)) for i in range(len(selected_csi))]

            for idx, csi_str in enumerate(selected_csi):
                tov_path = os.path.join(output_dir, model, b_str, csi_str, "tov.out")
                if os.path.exists(tov_path):
                    try:
                        data = np.loadtxt(tov_path)
                        if data.ndim < 2: continue
                        
                        # Columns: 0=e0, 1=Mass, 2=Radius
                        mess = data[:, 1]
                        radius = data[:, 2]
                        
                        # Label
                        try:
                            c_val = float(csi_str.replace("d","e"))
                            lbl = fr"$\xi={c_val:.1e}$"
                        except:
                            lbl = csi_str
                        
                        color = colors[idx]
                        line, = ax.plot(radius, mess, label=lbl, color=color, linewidth=1.2)
                        
                        # Max point
                        max_idx = np.argmax(mess)
                        ax.plot(radius[max_idx], mess[max_idx], marker='o', markersize=3, color=color)
                        
                        inset_curves.append((radius, mess, color))
                        
                    except Exception as e:
                        print(f"Skipping {tov_path}: {e}")

            ax.set_xlabel('Radius [km]')
            ax.set_xlim(9, 13.5)
            ax.set_ylim(bottom=1.0)
            ax.set_ylabel(r'Mass [$M_{\odot}$]')
            ax.set_title(f"{model} | $B={b_latex}$ G")
            ax.legend(fontsize='x-small', loc='center left', bbox_to_anchor=(1.0, 0.5))
            ax.grid(True, alpha=0.3, linestyle='--')
            
            # Call tight_layout before adding inset to avoid warnings
            plt.tight_layout()

            # --- Inset Zoom ---
            if inset_curves:
                ax_ins = inset_axes(ax, width="40%", height="40%", loc="lower left", 
                                    bbox_to_anchor=(0.1, 0.1, 1, 1), 
                                    bbox_transform=ax.transAxes)
                
                max_mass_global = 0
                max_rad_at_max = 0
                
                for r, m, c in inset_curves:
                    ax_ins.plot(r, m, color=c, linewidth=1.0)
                    m_max = np.max(m)
                    if m_max > max_mass_global:
                        max_mass_global = m_max
                        max_rad = r[np.argmax(m)]
                        max_rad_at_max = max_rad

                if max_mass_global > 0:
                    # Zoom window size
                    ax_ins.set_xlim(max_rad_at_max - 0.8, max_rad_at_max + 0.8)
                    ax_ins.set_ylim(max_mass_global - 0.05, max_mass_global + 0.005)
                
                ax_ins.tick_params(labelsize=6)
                mark_inset(ax, ax_ins, loc1=2, loc2=4, fc="none", ec="0.5", linewidth=0.5)
            
            fname = f"MR_{model}_{b_str}.pdf"
            plt.savefig(os.path.join(results_dir, fname), bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {fname}")

def plot_derivatives():
    """
    Plots the derivatives of M_max and R_at_max with respect to log10(csi).
    """
    results_dir = "results"
    csv_file = os.path.join(results_dir, "results_aggregate.csv")
    if not os.path.exists(csv_file):
        print(f"File not found: {csv_file}")
        return

    try:
        q = pl.scan_csv(csv_file)
        df = q.with_columns([
            pl.col("b").cast(pl.Utf8).str.replace("d", "e", literal=True).cast(pl.Float64, strict=False).alias("b_val"),
            pl.col("csi").cast(pl.Utf8).str.replace("d", "e", literal=True).cast(pl.Float64, strict=False).alias("csi_val"),
            pl.col("mmax").cast(pl.Float64),
            pl.col("rmax").cast(pl.Float64)
        ]).collect()
    except Exception as e:
        print(f"Error processing CSV: {e}")
        return

    if df.height == 0: return

    models = df["model"].unique().to_list()
    
    for model in models:
        mod_df = df.filter(pl.col("model") == model)
        b_vals = mod_df["b_val"].unique().sort().to_list()
        
        # Prepare figures
        fig_m, ax_m = plt.subplots(figsize=(6, 4))
        fig_r, ax_r = plt.subplots(figsize=(6, 4))
        
        has_data = False

        for b in b_vals:
            # Sort by csi
            subset = mod_df.filter(pl.col("b_val") == b).sort("csi_val")
            if len(subset) < 2: continue
            
            has_data = True
            csi = subset["csi_val"].to_numpy()
            mmax = subset["mmax"].to_numpy()
            rmax = subset["rmax"].to_numpy()
            
            # Filter csi > 0
            valid = csi > 0
            if not np.any(valid): continue
            
            csi = csi[valid]
            mmax = mmax[valid]
            rmax = rmax[valid]
            
            if len(csi) < 2: continue

            log_csi = np.log10(csi)
            
            # Derivatives
            d_log_csi = np.diff(log_csi)
            d_mmax = np.diff(mmax)
            d_rmax = np.diff(rmax)
            
            valid_diff = d_log_csi != 0
            if not np.any(valid_diff): continue

            deriv_m = d_mmax[valid_diff] / d_log_csi[valid_diff]
            deriv_r = d_rmax[valid_diff] / d_log_csi[valid_diff]
            
            # Midpoints
            x_mid = (log_csi[:-1] + log_csi[1:]) / 2
            x_plot = x_mid[valid_diff]

            # Labeling
            exp = int(np.log10(b)) if b > 0 else 0
            label_str = fr"$B=10^{{{exp}}}$ G"
            
            ax_m.plot(x_plot, deriv_m, linestyle='-', marker='o', markersize=3, alpha=0.7, label=label_str)
            ax_r.plot(x_plot, deriv_r, linestyle='-', marker='o', markersize=3, alpha=0.7, label=label_str)

        if has_data:
            # Save Derivative M plot
            ax_m.set_title(f"Model: {model} - d(Mmax)/d(log csi)")
            ax_m.set_xlabel(r"$\log_{10}(\xi)$")
            ax_m.set_ylabel(r"$d(M_{max})/d(\log_{10}(\xi))$")
            ax_m.legend(fontsize='x-small', title="Magnetic Field")
            ax_m.grid(True, alpha=0.3)
            
            out_root_m = f"deriv_mmax_vs_logcsi_{model}"
            plt.figure(fig_m.number)
            plt.tight_layout()
            plt.savefig(os.path.join(results_dir, f'{out_root_m}.pdf'))
            # plt.savefig(os.path.join(results_dir, f'{out_root_m}.png'), dpi=300)
            print(f"Saved: {out_root_m}")

            # Save Derivative R plot
            ax_r.set_title(f"Model: {model} - d(Rmax)/d(log csi)")
            ax_r.set_xlabel(r"$\log_{10}(\xi)$")
            ax_r.set_ylabel(r"$d(R_{at\_max})/d(\log_{10}(\xi))$")
            ax_r.legend(fontsize='x-small', title="Magnetic Field")
            ax_r.grid(True, alpha=0.3)
            
            out_root_r = f"deriv_rmax_vs_logcsi_{model}"
            plt.figure(fig_r.number)
            plt.tight_layout()
            plt.savefig(os.path.join(results_dir, f'{out_root_r}.pdf'))
            # plt.savefig(os.path.join(results_dir, f'{out_root_r}.png'), dpi=300)
            print(f"Saved: {out_root_r}")

        plt.close(fig_m)
        plt.close(fig_r)

if __name__ == "__main__":
    plot_m_vs_csi()
    plot_m_vs_rmax()
    plot_derivatives()
    plot_mr_curves_with_inset()
