import numpy as np
import nibabel as nib
import scipy.io

def save_seedmap_dscalar(
    seedmap: np.ndarray,
    template_cifti: nib.Cifti2Image,
    output_path: str,
    labels=None
):
    """
    seedmap shape: (n_grayordinates, n_maps)
    """

    if seedmap.ndim != 2:
        raise ValueError("Expected seedmap shape (grayordinates, n_maps)")

    n_gray, n_maps = seedmap.shape

    brain_models_axis = template_cifti.header.get_axis(1)

    if n_gray != brain_models_axis.size:
        raise ValueError(
            f"Grayordinate mismatch: {n_gray} vs {brain_models_axis.size}"
        )

    if labels is None:
        labels = [f"Seed_{i+1}" for i in range(n_maps)]

    if len(labels) != n_maps:
        raise ValueError("Number of labels must match number of maps")

    # CIFTI expects (maps, grayordinates)
    data = seedmap.T.astype(np.float32)

    scalar_axis = nib.cifti2.ScalarAxis(labels)

    header = nib.cifti2.Cifti2Header.from_axes(
        (scalar_axis, brain_models_axis)
    )

    img = nib.Cifti2Image(data, header)
    nib.save(img, output_path)
    print(f"Saved seedmap to {output_path}")

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

if __name__ == "__main__":

    raw_labels = [
        'DMN', 'Vis', 'FP', '',
        'DAN', '', 'VAN', 'Sal',
        'CO', 'SMd', 'SMl', 'Aud',
        'Tpole', 'MTL', 'PMN', 'PON', '', 'SCAN'
    ]

    labels = [
        name if name.strip() else f"Seed_{i+1}"
        for i, name in enumerate(raw_labels)
    ]

    #load seedmap from .mat file
    seedmap = load_seed_matrix('/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/seed_map_generation/all_data/python_seed_map_refractor_re/seedmaps_all_data_dtsereies_n_561_all_networks_Zscored.mat')
    template_cifti = nib.load('/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/template_matching_python/support_files/dscalar-template_data-mint_networks-15_colors-power_cifti-2.dscalar.nii')
    save_seedmap_dscalar(seedmap, template_cifti, '/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/seed_map_generation/all_data/python_seed_map_refractor_re/seedmaps_all_data_dtsereies_n_561_all_networks_Zscored.dscalar.nii', labels=labels)