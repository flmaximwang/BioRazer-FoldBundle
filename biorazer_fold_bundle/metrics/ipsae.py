import numpy as np


def calc_d0(L, pair_type="protein"):
    """
    Calculate the normalization parameter d0 for TM-score calculation.

    Parameters
    ----------
    L: int
        Length of the protein or nucleic acid sequence to align
    pair_type: str
        'protein' or 'nucleic_acid'
    """

    if pair_type == "protein":
        min_value = 1.0
    elif pair_type == "nucleic_acid":
        min_value = 2.0
    else:
        raise ValueError("pair_type must be 'protein' or 'nucleic_acid'")

    # When L < 27, d0 will be less than 1.0, which is not suitable for TM-score
    if L < 27:
        return min_value

    d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
    return max(min_value, d0)


def tm_score(x, d0):
    """
    Calculate TM-score based on
    """
    return 1.0 / (1.0 + (x / d0) ** 2)


def calc_ipsae(
    pae,
    atomarray,
    chain1,
    chain2,
    pae_cutoff=10.0,
    dist_cutoff=10.0,
    pair_type="protein",
):
    """
    Calculate ipSAE score for specified chain pairs

    Parameters
    ----------
    pae: numpy.ndarray, shape=(N, N)
        Pairwise predicted aligned error (PAE) matrix
    atomarray: biotite.structure.AtomArray
        AtomArray object containing atom coordinates and chain information corresponding to the PAE
    chain1, chain2: str
        target chains, like ('A', 'B')
    pae_cutoff, dist_cutoff: float
        cutoff for PAE and distance in atomarray. This makes sure only well predicted residues with close enough residues are analyzed
    pair_type: str
        'protein' / 'nucleic_acid'

    Return
    ------
    ipSAE: float
        The ipSAE score for the specified chain pair
    """
    pae = np.asarray(pae, dtype=np.float32)

    coords = atomarray.coord
    chains = atomarray.chain_id
    resnums = atomarray.res_id
    mask_ca = (
        atomarray.atom_name == "CA"
    )  # Because ipSAE is calculated based on residues
    coords = coords[mask_ca]
    chains = chains[mask_ca]
    resnums = resnums[mask_ca]
    idx1 = np.where(chains == chain1)[0]
    idx2 = np.where(chains == chain2)[0]
    if len(idx1) == 0 or len(idx2) == 0:
        raise ValueError(f"Chain {chain1} or {chain2} not found in the atom array.")

    dist_matrix = np.linalg.norm(
        coords[idx1][:, None, :] - coords[idx2][None, :, :], axis=-1
    )
    pae_sub = pae[np.ix_(idx1, idx2)]
    valid_mask = (pae_sub < pae_cutoff) & (
        dist_matrix < dist_cutoff
    )  # Considering only small pae and close residues

    ipsae_scores = np.zeros(len(idx1))

    if not np.any(valid_mask):
        return 0.0

    L = len(idx1) + len(idx2)
    d0 = calc_d0(L, pair_type)
    for i in range(len(idx1)):  # For every valid residue
        valid = valid_mask[i]  # Get residues in chain2 to align
        if np.any(valid):
            tm_score_tmp = tm_score(pae_sub[i][valid], d0)
            ipsae_scores[i] = np.mean(tm_score_tmp)
        else:
            ipsae_scores[i] = 0.0
    return max(ipsae_scores)  # Ensure non-negative score
