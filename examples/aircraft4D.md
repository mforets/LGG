We consider an LTI system that represents a longitudinal aircraft dynamics. The matrix $A$ represents stability derivatives of a Boeing 747 aircraft cruising at an altitude of 40kft with speed 774 ft/sec.

Let $\dot{x} = Ax + B u$, with state $x = [u,\alpha,\dot{\theta}, \theta]^T \in \mathbb{R}^4$ comprised of deviations in aircraft speed, angle of attack, pitch-rate, and pitch angle respectively. The matrix $A$ is 
$$
A = \begin{pmatrix}
-0.0030 & 0.0390 & 0 & -0.3220 \\
-0.0650 & -0.3190 & 7.7400 & 0 \\
0.0200 & -0.1010 & -0.4290 & 0 \\
0 & 0 & 1 & 0
\end{pmatrix},
$$
and 
$$
u = \begin{pmatrix}
0.0100 \\ -0.1800 \\ -1.1600 \\ 0
\end{pmatrix}
$$
The input is $u \in [-13.3,13.3] \in \mathbb{R}$ (measured in degrees).


---
References:
1. A. Bryson. Control of Spacecraft and Aircraft. Princeton Univ. Press, 1994.
2. S. Kaynama, M. Oishi. Schur-Based decomposition for reachability analysis of Linear-Time-Invariant systems. Joint 48th IEEE Conf on Decision and Control, 2009.
