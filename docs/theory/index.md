# Theory

This page summarizes the theoretical ideas behind the Alchemical Transfer Method (ATM) as implemented in AToM-OpenMM. It follows the treatment in Prof. Gallicchio's group note, [The Science Behind AToM](https://www.compmolbiophysbc.org/atom-openmm/the-science-behind-atom), with formatting adapted for this documentation.

## Alchemical Free Energy Perturbation

Alchemical free energy methods connect two thermodynamic states with a parameterized potential energy function. Let \(U_0(x)\) be the reference potential and \(U_1(x)\) be the perturbed potential for coordinates \(x\). The perturbation energy is

\[
u(x) = U_1(x) - U_0(x).
\]

An alchemical potential \(U_\lambda(x)\) interpolates between the two endpoints as the progress variable \(\lambda\) changes from 0 to 1. Samples collected at multiple alchemical states are combined with a free energy estimator, such as UWHAM, to estimate the free energy difference.

For binding, the reference state can be viewed as the ligand interacting with solvent, while the perturbed state corresponds to the ligand interacting with the receptor binding site.

## Alchemical Transfer

ATM differs from conventional alchemical approaches because it changes coordinates rather than force-field parameters. Instead of gradually turning interaction parameters on and off, ATM evaluates the energy change associated with transferring a ligand from one position to another.

For ABFE, the coordinate transformation translates the ligand from a solvent position into the binding site. If \(T\) denotes that transfer operation, the perturbed potential can be written schematically as

\[
U_1(x) = U_0(Tx),
\]

and the perturbation energy becomes

\[
u(x) = U_0(Tx) - U_0(x).
\]

The same idea extends to RBFE. One ligand is placed in solution and another in the binding site, and the alchemical coordinate transformation swaps their positions. ABFE and RBFE therefore use the same underlying idea: the physical meaning changes from translation of one ligand to coordinate swapping of two molecular partners.

Because the transformation is applied to coordinates, ATM can use standard molecular topologies and standard OpenMM energy functions. This is one of the central advantages of the method: the molecular mechanics potential does not need custom alchemical pair potentials to represent the endpoint states.

## Perturbation Function

AToM-OpenMM expresses the alchemical potential as

\[
U_\lambda(x) = U_0(x) + W_\lambda(u(x)),
\]

where \(W_\lambda(u)\) is the alchemical perturbation function. It must satisfy the endpoint conditions

\[
W_0(u) = 0,
\qquad
W_1(u) = u.
\]

A simple choice is the linear perturbation

\[
W_\lambda(u) = \lambda u.
\]

AToM-OpenMM commonly uses a softplus perturbation function instead:

\[
W_\lambda(u) =
\frac{\lambda_2-\lambda_1}{\alpha}
\ln \left[1+\exp\left(-\alpha(u-u_0)\right)\right]
+ \lambda_2 u + w_0.
\]

The schedule keywords `LAMBDA1`, `LAMBDA2`, `ALPHA`, `U0`, and `W0COEFF` define these parameters at each alchemical state. The softplus form behaves linearly at low and high perturbation energies, but it can use different slopes in those regimes. The linear perturbation is recovered when \(\lambda_1=\lambda_2=\lambda\).

This extra flexibility helps build smoother alchemical paths and improves overlap between neighboring states.

## Soft-Core Perturbation Energy

Alchemical binding calculations can suffer from endpoint singularities when atoms overlap. In ATM, large perturbation energies appear when a coordinate transfer places a ligand into a region already occupied by receptor or solvent atoms. Conventional approaches often address this with soft-core pair potentials.

AToM instead applies a soft-core mapping to the perturbation energy itself. The raw perturbation energy \(u\) is replaced by a bounded value \(u_\mathrm{sc}(u)\) before it is passed to \(W_\lambda\). In AToM-OpenMM, the relevant keywords are:

| Keyword | Role |
| --- | --- |
| `UMAX` | Upper bound approached by the soft-core perturbation energy. |
| `UBCORE` | Energy threshold below which the raw perturbation energy is left unchanged. |
| `ACORE` | Shape parameter controlling the smooth transition into the capped regime. |

For \(u \le U_\mathrm{bcore}\), the mapping leaves \(u\) unchanged. For larger values, it smoothly compresses the high-energy tail toward \(U_\mathrm{max}\). This preserves ordinary perturbation energies while preventing rare atomic-overlap samples from dominating averages.

The important conceptual point is that the soft-core transformation is applied to the scalar perturbation energy, not to every pair interaction in the force field. This keeps the method compatible with standard OpenMM potentials and avoids the need to re-evaluate stored configurations with a different lambda-dependent force field.

## Two Legs in Explicit Solvent

In explicit solvent, the full transformation is split into two alchemical legs. The split is used because a direct path all the way from one endpoint to the other can create severe solvent-ligand overlap near the opposite endpoint.

The first leg starts from the reference endpoint and proceeds to an alchemical intermediate near \(\lambda=1/2\). For a linear path, the intermediate has the symmetric form

\[
U_{1/2}(x) = \frac{1}{2}U_0(x) + \frac{1}{2}U_1(x).
\]

In this state, the transferred ligand has partial interaction with both environments. This avoids the most severe endpoint overlap problem while preserving a well-defined thermodynamic connection.

The second leg starts from the opposite endpoint and proceeds back to the same intermediate. If \(u(x)=U_1(x)-U_0(x)\), then the second leg uses the opposite perturbation, \(-u(x)\). If \(\Delta G_1\) and \(\Delta G_2\) are the free energy differences from the two endpoints to the intermediate, the binding free energy is obtained from their difference:

\[
\Delta G_b = \Delta G_1 - \Delta G_2.
\]

AToM-OpenMM marks the two legs with the `DIRECTION` schedule parameter. The `INTERMEDIATE` schedule parameter identifies states at the shared alchemical midpoint.

The soft-core mapping does not change the physical endpoint potentials. At the intermediate, the mapping should also leave relevant perturbation energies unchanged; in practice, this is checked by confirming that intermediate-state perturbation energies stay below the soft-core threshold.

## Relative Binding Free Energies

The same two-leg construction applies to RBFE. The endpoints correspond to opposite ligand placements: one ligand in solution and the other in the binding site, then the reverse placement. The intermediate is a nonphysical but useful alchemical state in which both partners partially interact with both environments.

AToM computes the perturbation energy by comparing energies before and after the coordinate swap. Because each energy evaluation places the ligands in different positions, the two ligands are not interpreted as directly interacting with each other in the intermediate state.

This is why ATM can treat ABFE and RBFE with a unified coordinate-transfer framework. Translation, swapping, and the variable-displacement workflows used by current examples are different coordinate transformations built around the same free energy machinery.

## Practical Connections

The theory maps directly onto the workflow keywords:

| Concept | AToM-OpenMM keywords |
| --- | --- |
| Alchemical states | `LAMBDAS`, `DIRECTION`, `INTERMEDIATE` |
| Softplus perturbation | `LAMBDA1`, `LAMBDA2`, `ALPHA`, `U0`, `W0COEFF` |
| Soft-core perturbation energy | `UMAX`, `ACORE`, `UBCORE` |
| Coordinate transfer | `DISPLACEMENT`, `LIGOFFSET`, ligand atom selections |
| Binding-site frame and restraints | `RCPT_CM_ATOMS`, `RCPT_FRAME_ATOMS_*`, `CM_KF`, `CM_TOL` |

For the workflow-level meaning of these keywords, see [Configuration Files](../user-guide/configuration.md).

## Source and Further Reading

This page is a documentation-oriented summary of [The Science Behind AToM](https://www.compmolbiophysbc.org/atom-openmm/the-science-behind-atom). The original page includes additional figures and CHK1 animations illustrating ligand transfer and coordinate swapping.
