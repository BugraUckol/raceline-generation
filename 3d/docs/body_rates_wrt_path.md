# Rates of Angular Errors Between the Body Frame & Path Frame

## Notation
* $s$ : Curvilinear progress of the closest point on the path
* $F_{E}$ : Earth Frame (Fixed)
* $F_{P}$ : Path Frame
* $F_{B}$ : Body Frame
* $\bar\omega_{B/E}^{(B)} = [p, q, r]^T $ : Angular velocity of the $F_{B}$ wrt. $F_{E}$ resolved in $F_{B}$
* $\bar\omega_{P/E}^{(P)} = \dot{s}[\tau, 0, \kappa]^T $ : Angular velocity of the $F_{P}$ wrt. $F_{E}$ resolved in $F_{P}$
* $\Phi = [e_\phi, e_\theta, e_\psi]^T$    
    * $e_\phi$ : Roll error of the $F_{B}$ wrt. $F_{P}$
    * $e_\theta$ : Pitch error of the $F_{B}$ wrt. $F_{P}$
    * $e_\psi$ : Yaw error of the $F_{B}$ wrt. $F_{P}$
* $ \hat{C}_{BE}^{(B)}$ : Transformation matrix that transforms elements of the $F_{B}$ to $F_{E}$ observed in the $F_{B}$
* $ \hat{C}_{PE}^{(P)}$ : Transformation matrix that transforms elements of the $F_{P}$ to $F_{E}$ observed in the $F_{P}$

## A Definition of The Angular Velocity
Given a transformation matrix $\hat{C}_{BE}^{(B)}$, $\bar\omega_{B/E}^{(B)}$ has the following relation

$\dot{\hat{C}}_{BE}^{(B)} = \hat{C}_{BE}^{(B)}\tilde\omega_{B/E}^{(B)} $

where the $\tilde\omega_{B/E}^{(B)}$ is the skew symmetric version of the $\bar\omega_{B/E}^{(B)}$

## Relative Angular Velocity Between Two Frames

$\vec\omega_{B/P} = \vec\omega_{B/E} - \vec\omega_{P/E}$

This vector equation can be resolved in any frame. However, input of the system is easier to express on the body frame, thus, $F_B$ is the natural candidate. Since $\bar\omega_{P/E}^{(B)}$ do not exist directly it must be found using the relation below (Equation 3.2.9 in ME502 lecture notes)

$\bar\omega_{P/E}^{(B)} = \hat{C}_{BP}^{(B),T} \bar\omega_{P/E}^{(P)}$

Then the formula takes the following form

$\bar\omega_{B/P}^{(B)} = \bar\omega_{B/E}^{(B)} -\hat{C}_{BP}^{(B),T} \bar\omega_{P/E}^{(P)}$

```
Note: Transformation dyadic between two frames will have the same matrix representation when observed in either one of the frames
```

## Euler Angle Rates Having The Relative Angular Velocity
Angular velocity, when projected and subtracted in the intermediate frames, correspond to the Euler Angles' rates. This is due to the total relative angular velocity is the sum of the all three of the relative angular velocities in between.

This operation can be written as a simple matrix operation, that is well known.

$\dot{{\Phi}} = \hat H \bar\omega_{B/P}^{(B)}$

## Computer Programming
All the functions that are necessary are located at the Rigid Body Dynamics repository.
