import pandas
import seaborn
import allel


def do_pca(
    genotypes,
    ld_pruning_params=None,
    n_components=10,
    scaler=None,
    use_randomized_pca_optimization=False,
):
    if ld_pruning_params is None:
        ld_pruning_params = {}

    num_variants_before_ld_pruning = genotypes.num_variants
    genotypes = genotypes.ld_prune(**ld_pruning_params)
    num_variants_after_ld_pruning = genotypes.num_variants
    genotypes = genotypes.to_genotype_012_matrix()

    if use_randomized_pca_optimization:
        pca_funct = allel.randomized_pca
    else:
        pca_funct = allel.pca
    projections, model = pca_funct(genotypes, n_components=n_components, scaler=scaler)

    explained_variances = model.explained_variance_ratio_
    princom_names = [f"princomp_{idx + 1}" for idx in range(explained_variances.size)]
    explained_variances = pandas.Series(explained_variances, index=princom_names)
    explained_variances = explained_variances * 100.0

    projections = pandas.DataFrame(projections, columns=princom_names)

    return {
        "projections": projections,
        "explained_variances": explained_variances,
        "num_variants_after_ld_pruning": num_variants_after_ld_pruning,
        "num_variants_before_ld_pruning": num_variants_before_ld_pruning,
    }


def plot_pca_result(
    pca_result, axes, axis1="princomp_1", axis2="princomp_2", classification=None
):
    seaborn.scatterplot(
        data=pca_result["projections"], x=axis1, y=axis2, hue=classification, ax=axes
    )
    explained_variances = pca_result["explained_variances"]
    var1 = explained_variances[axis1]
    var2 = explained_variances[axis2]
    axes.set_xlabel(f"{axis1} ({var1:.1f}%)")
    axes.set_ylabel(f"{axis2} ({var2:.1f}%)")
