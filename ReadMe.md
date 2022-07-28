# Compare different nonlinear filters

### Available filters : 
 - particle filter
 - extended Kalman filter
 - central-difference Kalman filter



# Central-difference Kalman filter

## Reference:

  - King2021 - Struktur- und Parameteridentifikation.
  - Noergaard2000 - New developments in state estimation for nonlinear systems. https://www.sciencedirect.com/science/article/pii/S0005109800000893

 

## System desciption: 

$$\begin{align}
    x(k+1)   &= f(x(k),u(k)) + w(k) \\ 
            y(k) &= Cx(k) + v(k) 
\end{align}$$

## Assumptions: 
- nonlinear ODE with linear output function
- process and measurement noise are zero-mean, additive and Gaussian 
   - $w_k \sim \mathcal{N}(0,\,Q)$ , $v_k \sim \mathcal{N}(0,\,R)$

## Algorithm:

### 1. Initialization
 - $\hat{x}(0|0), P_{xx}(0|0), Q, R, h$
 - h is the selected interval length of the CD-approaximation (set $h = \sqrt{3}$ for Gaussian distribution)

### 2. Time-Update / Prediction
Definítion of sigma points for x: 

$$\begin{align*}
\mathcal{X}_{0}(k+1|k) &=  f\left(\hat{x}(k|k) ,u_k \right)\\
\mathcal{X}_{plus,i}(k+1|k) &=   f\left(\hat{x}(k|k) + h\cdot \left(\sqrt{P_{xx}(k|k)}\right)_i ,u_k \right),  i = 1,...,n_x\\
\mathcal{X}_{minus,i}(k+1|k) &= f\left(\hat{x}(k|k) - h\cdot \left(\sqrt{P_{xx}(k|k)}\right)_i ,u_k \right),  i = 1,...,n_x
\end{align*}$$

- $\left(\sqrt{P}\right)_i$ is the i-th coloum of Cholesky decomposition of $P$, i.e. i-th column of $S$ with $P = SS^T$
   - Notice: the result of MATLAB chol( ) command is the transpose

Predicted state: 

$$\begin{align*}
    \hat{x}(k+1|k)   &= \frac{h^2-n_x}{h^2} \mathcal{X}_{0}(k+1|k) \\ 
             &+ \frac{1}{2h^2} \sum^{n_x}_{i=1} \left[\mathcal{X}_{plus,i}(k+1|k)\\
             + \mathcal{X}_{minus,i}(k+1|k)  \right]
\end{align*}$$

- The influence of $w_k$ is concealed, since it is indeed linear (additive after nonlinear time propogation).


Predicted covariance matrix of the state: 

$$\begin{align*}
    P_{xx}(k+1|k)   =& E\left( x(k+1) -\hat{x}(k+1|k)  \right)\left( x(k+1) -\hat{x}(k+1|k)  \right)^T \\
    =& \frac{1}{4h^2} \sum^{n_x}_{i=1} \left[\mathcal{X}_{plus,i}(k+1|k)\\
     - \mathcal{X}_{minus,i}(k+1|k)  \right]^2\\
    &+ \frac{h^2-1}{4h^2} \sum^{n_x}_{i=1} \left[\mathcal{X}_{plus,i}(k+1|k)\\
     + \mathcal{X}_{minus,i}(k+1|k) - 2\mathcal{X}_{0}(k+1|k)  \right]^2 \\
    & + Q(k+1)         
\end{align*}$$

```sh
# Replace for-loop by matrix multiplication
      tmpVar1 = sigmaXi_kplus1_k_plus - sigmaXi_kplus1_k_minus;
      tmpVar2 = sigmaXi_kplus1_k_plus + sigmaXi_kplus1_k_minus - 2*sigmaX0_kplus1_k;
      tmpVar3 = [sqrt(1/(4*h^2))*tmpVar1  sqrt((h^2-1)/(4*h^4))*tmpVar2];
      Pxx_kplus1_k = tmpVar3*transpose(tmpVar3) + Q; 
```




### 3.1 Measurement-Update / Correction  (linear output function)
Predicted output:

$$\begin{align*}
    \hat{y}(k+1|k)   &= C\hat{x}(k+1|k)
\end{align*}$$
 - derivation see King2021-SPI


Predicted cross covariance matrix:

$$\begin{align*}
    P_{xy}(k+1|k)   =& E\left( x(k+1) -\hat{x}(k+1|k)  \right)\left( y(k+1) -\hat{y}(k+1|k)  \right)^T \\
    =&     E\left( x(k+1) -\hat{x}(k+1|k)  \right)
     \left(C (x(k+1) -\hat{x}(k+1|k)) + v(k+1)  \right)^T \\
    =&     E\left( x(k+1) -\hat{x}(k+1|k)  \right)\left( x(k+1) -\hat{x}(k+1|k)  \right)^T C^T \\
    & +  E\left( x(k+1) -\hat{x}(k+1|k)  \right)v(k+1)^T  \\
    =& E\left( x(k+1) -\hat{x}(k+1|k)  \right)\left( x(k+1) -\hat{x}(k+1|k)  \right)^T C^T \\
    =& P_{xx}(k+1|k)C^T
\end{align*}$$

Predicted covariance matrix of the output:

$$\begin{align*}
    P_{yy}(k+1|k)   =& E\left( y(k+1) -\hat{y}(k+1|k)  \right)\left( y(k+1) -\hat{y}(k+1|k)  \right)^T \\
    =& CP_{xx}(k+1|k)C^T + R(k+1)
\end{align*}$$
 - The derivation is trivial :-)


Correction gain: 

$$\begin{align*}
    K(k+1)   = P_{xy}(k+1|k) P^{-1}_{yy}(k+1|k)
\end{align*}$$


Corrected esimation: 

$$\begin{align*}
    \hat{x}(k+1|k+1)   &= \hat{x}(k+1|k) + K(k+1)\left( y(k+1) - \hat{y}(k+1|k) \right) \\
    P_{xx}(k+1|k+1) &= P_{xx}(k+1|k) - K(k+1)P_{yy}(k+1|k)K(k+1)^T
\end{align*}$$



### 3.2 Measurement-Update / Correction  (nonlinear output function)
Nonlinear output function: $y(k) = g(x(k)) + v(k)$

Definítion of sigma points for y: 

$$\begin{align*}
\mathcal{Y}_{0}(k+1|k) &=   
             g\left(\hat{x}(k+1|k) \right)\\
\mathcal{Y}_{plus,i}(k+1|k) &=   
             g\left(\hat{x}(k+1|k) + h\cdot \left(\sqrt{P_{xx}(k+1|k)}\right)_i \right),  i = 1,...,n_y\\
\mathcal{Y}_{minus,i}(k+1|k) &=
             g\left(\hat{x}(k+1|k) - h\cdot \left(\sqrt{P_{xx}(k+1|k)}\right)_i  \right),  i = 1,...,n_y
\end{align*}$$

- $\left(\sqrt{P}\right)_i$ is the i-th coloum of Cholesky decomposition of $P$, i.e. i-th column of $S$ with $P = SS^T$
   - Notice: the result of MATLAB chol( ) command is the transpose


Predicted output:

$$\begin{align*}
    \hat{y}(k+1|k)   &= \frac{h^2-n_x}{h^2} \mathcal{Y}_{0}(k+1|k) \\ 
             &+ \frac{1}{2h^2} \sum^{n_x}_{i=1} \left[
             \mathcal{X}_{plus,i}(k+1|k)\\
             + 
             \mathcal{X}_{minus,i}(k+1|k)  \right]
\end{align*}$$

- The influence of $w_k$ is concealed, since it is indeed linear (additive after nonlinear time propogation).


Predicted cross covariance matrix:

$$\begin{align*}
    P_{xy}(k+1|k)   = \frac{1}{2h} \sqrt{P_{xx}(k+1|k)} 
    \left(\mathcal{Y}_{plus,1:n_y}(k+1|k)) - \mathcal{Y}_{minus,1:n_y}(k+1|k) \right)^T    
\end{align*}$$

Predicted covariance matrix of the output:

$$\begin{align*}
    P_{yy}(k+1|k)   =& E\left( y(k+1) -\hat{y}(k+1|k)  \right)\left( y(k+1) -\hat{y}(k+1|k)  \right)^T \\
    =& \frac{1}{4h^2} \sum^{n_x}_{i=1} \left[
             \mathcal{Y}_{plus,i}(k+1|k)\\
             - 
             \mathcal{Y}_{minus,i}(k+1|k)  \right]^2
             \\
    &+ \frac{h^2-1}{4h^2} \sum^{n_x}_{i=1} \left[
             \mathcal{Y}_{plus,i}(k+1|k)\\
             + 
             \mathcal{Y}_{minus,i}(k+1|k) 
             - 2\mathcal{Y}_{0}(k+1|k)  \right]^2 \\
    & + R(k+1) 
\end{align*}$$

```sh
# Replace for-loop by matrix multiplication
      tmpVar1 = sigmaYi_kplus1_k_plus - sigmaYi_kplus1_k_minus;
      tmpVar2 = sigmaYi_kplus1_k_plus + sigmaYi_kplus1_k_minus - 2*sigmaY0_kplus1_k;
      tmpVar3 = [sqrt(1/(4*h^2))*tmpVar1  sqrt((h^2-1)/(4*h^4))*tmpVar2];
      Pyy_kplus1_k = tmpVar3*transpose(tmpVar3) + Q; 
```

Correction gain: 

$$\begin{align*}
    K(k+1)   = P_{xy}(k+1|k) P^{-1}_{yy}(k+1|k)
\end{align*}$$


Corrected esimation: 

$$\begin{align*}
    \hat{x}(k+1|k+1)   &= \hat{x}(k+1|k) + K(k+1)\left( y(k+1) - \hat{y}(k+1|k) \right) \\
    P_{xx}(k+1|k+1) &= P_{xx}(k+1|k) - K(k+1)P_{yy}(k+1|k)K(k+1)^T
\end{align*}$$

