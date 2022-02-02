import numpy
import pandas


def _center_matrix(matrix):
    "It centers the matrix"
    means = matrix.mean(axis=0)
    return matrix - means


def do_pca_from_genotypes(genotypes_012_with_snps_in_rows: pandas.DataFrame):

    "It does a Principal Component Analysis"

    # transform the genotype data into a 2-dimensional matrix where each cell
    # has the number of non-reference alleles per call

    matrix = genotypes_012_with_snps_in_rows.T

    n_rows, n_cols = matrix.shape
    if n_cols < n_rows:
        # This restriction is in the matplotlib implementation, but I don't
        # know the reason
        msg = "The implementation requires more SNPs than samples"
        raise RuntimeError(msg)

    # Implementation based on the matplotlib PCA class
    cen_matrix = _center_matrix(matrix)
    # The following line should be added from a example to get the correct
    # variances
    # cen_scaled_matrix = cen_matrix / math.sqrt(n_rows - 1)
    cen_scaled_matrix = cen_matrix

    singular_vals, princomps = numpy.linalg.svd(cen_scaled_matrix, full_matrices=False)[
        1:
    ]

    samples = genotypes_012_with_snps_in_rows.columns
    princom_names = [f"princomp_{idx + 1}" for idx in range(princomps.shape[0])]
    princomps = pandas.DataFrame(princomps, index=samples)

    eig_vals = singular_vals ** 2

    eig_vals = pandas.Series(eig_vals, index=princom_names)
    var_percentages = eig_vals / eig_vals.sum() * 100.0
    projections = numpy.dot(princomps, cen_matrix.T).T
    projections = pandas.DataFrame(projections, index=samples, columns=princom_names)

    return {
        "eigenvalues": eig_vals,
        "projections": projections,
        "var_percentages": var_percentages,
        "principal_components": princomps,
    }
