# Rates of Angular Errors Between the Body Frame & Path Frame
## System Description
Consider a rigid object whose velocity vector is inline with its x-axis and has a constant magnitude. Point that is closest to it on the path has an associated Frenet Frame about which the equations are to be written.
## Notation
### Frames
* $s$ : Curvilinear progress of the closest point on the path
* $F_{E}$ : Earth Frame (Fixed)
* $F_{P}$ : Path Frame
* $F_{B}$ : Body Frame
### Positions
* $\vec R_{B/E}$ : Position vector of the B wrt. origin of E
* $\vec R_{P/E}$ : Position vector of the P wrt. origin of E
* $\vec R_{B/P}$ : Position vector of the B wrt. origin of P
### Orientation $\Phi = [e_\phi, e_\theta, e_\psi]^T$    
* $e_\phi$ : Roll error of the $F_{B}$ wrt. $F_{P}$
* $e_\theta$ : Pitch error of the $F_{B}$ wrt. $F_{P}$
* $e_\psi$ : Yaw error of the $F_{B}$ wrt. $F_{P}$
* $ \hat{C}_{BE}^{(B)}$ : Transformation matrix that transforms elements of the $F_{B}$ to $F_{E}$ observed in the $F_{B}$
* $ \hat{C}_{PE}^{(P)}$ : Transformation matrix that transforms elements of the $F_{P}$ to $F_{E}$ observed in the $F_{P}$
### Angular Velocities
* $\bar\omega_{B/E}^{(B)} = \begin{bmatrix}
                p \\
                q \\
                r
                \end{bmatrix} $ : Angular velocity of the $F_{B}$ wrt. $F_{E}$ resolved in $F_{B}$
* $\bar\omega_{P/E}^{(P)} = \dot{s}\begin{bmatrix}
                \tau\\
                0 \\
                \kappa
                \end{bmatrix} $ : Angular velocity of the $F_{P}$ wrt. $F_{E}$ resolved in $F_{P}$

## Kinematical Relations
### Resolution of Known Vectors
$\bar R_{B/P}^{(P)} = \begin{bmatrix}
                    0\\
                    e_n \\
                    e_b
                    \end{bmatrix}$

$\bar V^{(P)} = 
\begin{bmatrix}
  cos(e_{\theta})cos(e_{\psi})\\
  cos(e_{\theta})sin(e_{\psi}) \\
  -sin(e_{\theta})
\end{bmatrix}$

$D_E\vec R_{P/E}^{(P)} = \begin{bmatrix}
                    \dot s\\
                    0 \\
                    0
                    \end{bmatrix}$

### Loop closure

$\vec R_{B/E} = \vec R_{P/E} + \vec R_{B/P}$

Taking the derivative wrt. the fixed frame

$D_E\vec R_{B/E} = D_E\vec R_{P/E} + D_E\vec R_{B/P}$

Using Coriolis Transport Theorem

$\vec V = D_E\vec R_{P/E} + D_P\vec R_{B/P} + \bar\omega_{P/E} \times \vec R_{B/P}$

Resolving this equation in $F_{P}$

$V\begin{bmatrix}
  cos(e_{\theta})cos(e_{\psi})\\cos(e_{\theta})sin(e_{\psi})\\-sin(e_{\theta})\end{bmatrix} = \begin{bmatrix}0\\\dot e_n \\\dot e_b\end{bmatrix} + \begin{bmatrix}\dot s\\0 \\0\end{bmatrix} + \dot s\begin{bmatrix}- e_n \kappa\\-e_b \tau\\e_n \tau\end{bmatrix}$

Factoring $\dot s$ out

$V\begin{bmatrix}
  cos(e_{\theta})cos(e_{\psi})\\cos(e_{\theta})sin(e_{\psi})\\-sin(e_{\theta})\end{bmatrix} = \begin{bmatrix}0\\\dot e_n \\\dot e_b\end{bmatrix} + \dot s\begin{bmatrix} 1 - e_n \kappa\\-e_b \tau\\e_n \tau\end{bmatrix}$

Notice that $\dot s$ equation can be written separately

$\dot s = \frac{Vcos(e_{\theta})cos(e_{\psi})}{1 - e_n \kappa}$ (1)

Then the dynamics for lateral and vertical error can be written as follows

$\dot e_n = Vcos(e_{\theta})sin(e_{\psi}) + \dot s \tau e_b$ (2)

$\dot e_b = -Vsin(e_{\theta}) - \dot s \tau e_n$ (3)

These three equations will form the basis of the movement of the $F_{P}$ and the relative errors of the $F_{B}$ wrt. the $F_{P}$. By controlling the orientation of the velocity vector (x-axis of the $F_{B}$) one might control the lateral and vertical deviations. 