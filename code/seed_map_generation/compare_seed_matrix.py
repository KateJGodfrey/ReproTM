import scipy.io
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt
import seaborn as sns

def load_seed_matrix(filepath):
    """Load seed_matrix from a .mat file."""
    try:
        mat_data = scipy.io.loadmat(filepath)
        if 'seed_matrix' not in mat_data:
            raise KeyError(f"'seed_matrix' not found in {filepath}")
        return mat_data['seed_matrix']
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        raise

def compute_similarity_metrics(matrix1, matrix2):
    """Compute multiple similarity metrics between two matrices."""
    results = {}
    
    # Check if matrices have the same shape
    if matrix1.shape != matrix2.shape:
        print(f"ERROR: Matrices must have the same shape!")
        print(f"Matrix 1: {matrix1.shape}, Matrix 2: {matrix2.shape}")
        return None
    
    # Flatten matrices
    flat1 = matrix1.flatten()
    flat2 = matrix2.flatten()
    
    # 1. Pearson Correlation (measures linear relationship)
    pearson_r, pearson_p = pearsonr(flat1, flat2)
    results['pearson_r'] = pearson_r
    results['pearson_p'] = pearson_p
    
    # 2. Spearman Correlation (measures monotonic relationship, rank-based)
    spearman_r, spearman_p = spearmanr(flat1, flat2)
    results['spearman_r'] = spearman_r
    results['spearman_p'] = spearman_p
    
    # 3. Cosine Similarity (measures orientation similarity)
    cosine_sim = np.dot(flat1, flat2) / (np.linalg.norm(flat1) * np.linalg.norm(flat2))
    results['cosine_similarity'] = cosine_sim
    
    # 4. Mean Squared Error (MSE) - lower is better
    mse = mean_squared_error(flat1, flat2)
    results['mse'] = mse
    
    # 5. Root Mean Squared Error (RMSE) - lower is better
    rmse = np.sqrt(mse)
    results['rmse'] = rmse
    
    # 6. Mean Absolute Error (MAE) - lower is better
    mae = mean_absolute_error(flat1, flat2)
    results['mae'] = mae
    
    # 7. Normalized RMSE (as percentage of data range)
    data_range = max(flat1.max(), flat2.max()) - min(flat1.min(), flat2.min())
    nrmse = (rmse / data_range) * 100 if data_range != 0 else 0
    results['nrmse_percent'] = nrmse
    
    # 8. R-squared (coefficient of determination)
    ss_res = np.sum((flat1 - flat2) ** 2)
    ss_tot = np.sum((flat1 - np.mean(flat1)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
    results['r_squared'] = r_squared
    
    # 9. Maximum absolute difference
    max_diff = np.max(np.abs(matrix1 - matrix2))
    results['max_abs_diff'] = max_diff
    
    # 10. Frobenius norm of difference (matrix-specific)
    frobenius_norm = np.linalg.norm(matrix1 - matrix2, 'fro')
    results['frobenius_norm'] = frobenius_norm
    
    # 11. Percentage of elements within tolerance
    tolerance = 0.01 * data_range  # 1% of data range
    within_tolerance = np.sum(np.abs(matrix1 - matrix2) < tolerance) / matrix1.size * 100
    results['within_1pct_tolerance'] = within_tolerance
    
    return results

def visualize_comparison(matrix1, matrix2, save_path=None):
    """Create visualization comparing the two matrices."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot matrix 1
    im1 = axes[0, 0].imshow(matrix1, cmap='viridis', aspect='auto')
    axes[0, 0].set_title('Matrix 1')
    plt.colorbar(im1, ax=axes[0, 0])
    
    # Plot matrix 2
    im2 = axes[0, 1].imshow(matrix2, cmap='viridis', aspect='auto')
    axes[0, 1].set_title('Matrix 2')
    plt.colorbar(im2, ax=axes[0, 1])
    
    # Plot difference
    diff = matrix1 - matrix2
    im3 = axes[1, 0].imshow(diff, cmap='RdBu_r', aspect='auto')
    axes[1, 0].set_title('Difference (Matrix1 - Matrix2)')
    plt.colorbar(im3, ax=axes[1, 0])
    
    # Scatter plot
    axes[1, 1].scatter(matrix1.flatten(), matrix2.flatten(), alpha=0.5, s=1)
    axes[1, 1].plot([matrix1.min(), matrix1.max()], 
                     [matrix1.min(), matrix1.max()], 
                     'r--', label='Perfect match')
    axes[1, 1].set_xlabel('Matrix 1 values')
    axes[1, 1].set_ylabel('Matrix 2 values')
    axes[1, 1].set_title('Value-by-value comparison')
    axes[1, 1].legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\nVisualization saved to: {save_path}")
    else:
        plt.show()

def print_results(results):
    """Print similarity metrics in a formatted way."""
    print("\n" + "="*70)
    print("SIMILARITY METRICS")
    print("="*70)
    
    print("\n1. CORRELATION METRICS (higher is better, range: -1 to 1)")
    print(f"   Pearson Correlation:    {results['pearson_r']:>8.4f}  (p={results['pearson_p']:.4e})")
    print(f"   Spearman Correlation:   {results['spearman_r']:>8.4f}  (p={results['spearman_p']:.4e})")
    print(f"   Cosine Similarity:      {results['cosine_similarity']:>8.4f}")
    print(f"   R-squared:              {results['r_squared']:>8.4f}")
    
    print("\n2. ERROR METRICS (lower is better)")
    print(f"   Mean Squared Error:     {results['mse']:>8.4f}")
    print(f"   Root Mean Squared Error:{results['rmse']:>8.4f}")
    print(f"   Mean Absolute Error:    {results['mae']:>8.4f}")
    print(f"   Normalized RMSE:        {results['nrmse_percent']:>8.2f}%")
    print(f"   Max Absolute Diff:      {results['max_abs_diff']:>8.4f}")
    print(f"   Frobenius Norm:         {results['frobenius_norm']:>8.4f}")
    
    print("\n3. SIMILARITY ASSESSMENT")
    print(f"   % within 1% tolerance:  {results['within_1pct_tolerance']:>8.2f}%")
    
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    
    # Overall assessment
    pearson = results['pearson_r']
    nrmse = results['nrmse_percent']
    within_tol = results['within_1pct_tolerance']
    
    if pearson > 0.99 and nrmse < 1 and within_tol > 99:
        print("✓ EXCELLENT: Matrices are nearly identical")
    elif pearson > 0.95 and nrmse < 5 and within_tol > 90:
        print("✓ VERY GOOD: Matrices are very similar")
    elif pearson > 0.90 and nrmse < 10:
        print("~ GOOD: Matrices are similar with minor differences")
    elif pearson > 0.80:
        print("~ MODERATE: Matrices show moderate similarity")
    elif pearson > 0.60:
        print("⚠ FAIR: Matrices have notable differences")
    else:
        print("✗ POOR: Matrices are substantially different")
    
    print("="*70)

def main():
    # File paths - modify these to your actual file paths
    file1 = '/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/seed_map_generation/all_data/python_seed_map_refractor_re/seedmaps_all_data_dtsereies_n_561_all_networks_Zscored.mat'
    file2 = '/projects/standard/rando149/shared/projects/ABCD_template_maker/seedmaps_from_template2.0/ABCD_scan_seedmap/seedmaps_subs_withsmoothed_dtseries_n141_all_networksZscored.mat'
    
    print("Loading matrices...")
    matrix1 = load_seed_matrix(file1)
    matrix2 = load_seed_matrix(file2)
    
    print(f"\nMatrix 1 shape: {matrix1.shape}")
    print(f"Matrix 2 shape: {matrix2.shape}")

    # Check and fix NaN values in matrix2 column 16
    if matrix2.shape[1] > 16:  # Column 16 is index 15
        col_16 = matrix2[:, 16]
        if np.any(np.isnan(col_16)):
            print("Warning: NaN values found in matrix2 column 16. Replacing with 0...")
            matrix2[:, 16] = np.nan_to_num(col_16, nan=0.0)
    
    print(f"\nMatrix 1 statistics:")
    print(f"  Mean: {np.mean(matrix1):.4f}")
    print(f"  Std: {np.std(matrix1):.4f}")
    print(f"  Min: {np.min(matrix1):.4f}")
    print(f"  Max: {np.max(matrix1):.4f}")
    
    print(f"\nMatrix 2 statistics:")
    print(f"  Mean: {np.mean(matrix2):.4f}")
    print(f"  Std: {np.std(matrix2):.4f}")
    print(f"  Min: {np.min(matrix2):.4f}")
    print(f"  Max: {np.max(matrix2):.4f}")
    
    print("\nComparing matrices...")
    # Compute all similarity metrics
    results = compute_similarity_metrics(matrix1, matrix2)
    
    if results is None:
        return
    
    # Print results
    print_results(results)
    
    # Create visualization
    print("\nGenerating comparison visualization...")
    visualize_comparison(matrix1, matrix2, save_path='/scratch.global/pandh015/precision_maps_via_template_matching/sanju_vs_prod_TM/matrix_comparison_fulldata_vs_sanju_map.png')

if __name__ == "__main__":
    main()