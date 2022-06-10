import numpy as np

# ----------------------------------------------------------------------------------------
def local_branching(tree, tau=1, tau0=1, infinite_root_branch=True):
    r"""Add local branching statistics (Neher et al. 2014) as tree node
    features to the ETE tree attribute.
    After execution, all nodes will have new features ``LBI``
    (local branching index) and ``LBR`` (local branching ratio, below Vs
    above the node)

    Args:
        tau: decay timescale for exponential filter
        tau0: effective branch length for branches with zero mutations
        infinite_root_branch: calculate assuming the root node has an infinite branch
    """
    # the fixed integral contribution for clonal cells indicated by abundance annotations
    clone_contribution = tau * (1 - np.exp(-tau0 / tau))


    # post-order traversal to populate downward integrals for each node
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf():
            node.LB_down = {
                node: node.abundance * clone_contribution
                if node.abundance > 1
                else 0
            }
        else:
            node.LB_down = {node: node.abundance * clone_contribution}
            for child in node.children:
                node.LB_down[child] = tau * (
                    1 - np.exp(-child.dist / tau)
                ) + np.exp(-child.dist / tau) * sum(child.LB_down.values())


    # pre-order traversal to populate upward integral for each node
    for node in tree.traverse(strategy="preorder"):
        if node.is_root():
            # integral corresponding to infinite branch above root node
            node.LB_up = tau if infinite_root_branch else 0
        else:
            node.LB_up = tau * (1 - np.exp(-node.dist / tau)) + np.exp(
                -node.dist / tau
            ) * (
                node.up.LB_up
                + sum(
                    node.up.LB_down[message]
                    for message in node.up.LB_down
                    if message != node
                )
            )


    # finally, compute LBI (LBR) as the sum (ratio) of downward and upward integrals at each node
    for node in tree.traverse():
        node_LB_down_total = sum(node.LB_down.values())
        node.LBI = node_LB_down_total + node.LB_up
        node.LBR = node_LB_down_total / node.LB_up if not node.is_root() else 0.
