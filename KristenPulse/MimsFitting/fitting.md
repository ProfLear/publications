# Fitting

This file contains the information on the types of fitting I have tried. The core of the fitting routine is two expressions that calculates splitting constants for hydrogens as a function of distance from the nanoparticle surface and angle from the applied magnetif field. These are for contact and dipolar coupling.

### Contact coupling

$A_{\textrm{contact}}(r) = \frac{2\mu_0 g_e \beta_e g_n \beta_n}{3h}(I_0e^{-\frac{r - r_{np}}{\ell}})^2$

In this equation $\mu_0$ is the permeability of free space, $g_e$ is the electron g-value, $\beta_e$ is the Bohr magneton, $g_n$ is the nuclear g-value, $\beta_n$ is the nuclear magneton, $h$ is Planck's constant, $I_0$ is the amplitude of the wavefunction at the surface of the nanoparticle, $r$ is the distance from the nanoparticle center, $r_{np}$ is the radius of the nanoparticle, and $\ell$ is the characteristic decay length for the wavefunction from the surface of the nanoparticle.

It is also true, however, that we might expect $I_0$ to depend on the radius of the nanoparticle.

For a spherical nanoparticle, we assume the electron density exists primarily outside the particle (tunneling into the surrounding matrix). The normalization integral in spherical coordinates is:

‬‭$\int_{R_{np}}^{\infty} |\psi(r)|^2 \cdot 4\pi r^2 \, dr = 1$‬

Given the above form of the wave function ‭‬‭‬‭‬‭‬‭‬, the integral becomes:

‭$4\pi I_0^2 \int_{R_{np}}^{\infty} e^{-2(r - R_{np})/\ell} \cdot r^2 \, dr = 1$

‬
By performing a change of variables (‭‬‭‬‭‬‬‭‬‭‬‬‭‬‭‬$x = r - r_{np}$) and solving the resulting integral, we find the relationship for ‭‬:

‭$I_0(r_{np}) = \sqrt{\frac{1}{2\pi \ell (R_{np}^2 + R_{np}\ell + \frac{1}{2}\ell^2)}}$

‬This shows that, as the nanoparticle size increases, we wave function spreads out and the amplitude at the surface decreases.

So, the total contact coupling becomes:

$A_{\textrm{contact}}(r) = \frac{2}{3}\frac{\mu_0 g_e \beta_e g_n \beta_n}{h}(‭\sqrt{\frac{1}{2\pi \ell (R_{np}^2 + R_{np}\ell + \frac{1}{2}\ell^2)}})(e^{-\frac{r - r_{np}}{\ell}})^2$

### Dipolar coupling

$A_\textrm{dipolar}(r) = \frac{1}{4\pi}\frac{\mu_0 g_e \beta_e g_n \beta_n}{h}(\frac{1}{r^3})(3 \cos^2\theta - 1)$

These produce a coupling constant that is in units of Hz.

You can see in the above that both equations share a common factor of $\frac{\mu_0 g_e \beta_e g_n \beta_n}{h}$

## General thoughts

There are a number of fitting routines that I want to try.  Each will look at coupling shell by shell away from the nanoparticle, using information on hydrogen distributions from molecular dynamics simulations. We will need to supply infomration about the nanoparticle size, the electron and nuclear g-factor and so on.

I also want to be able to check if Lorentzian or Gaussian lineshapes are better.

So, overall:

1. Dependence on coupling mechanism (dipolar, contact, or both)
2. Dependence on size
3. Dependence on lineshape
4. Dependence on broadening strain
5. ???

## Mean size nanoparticles.

For these calculations, the idea is to take the mean value that came from the fitting of the nanoparticle distribtutions. For this, the resulting lineshape of the individual features might need to be Gaussian? Will need to check.

### Dipolar only

No contact coupling.  So, this handles the case that there is no appreciable effects coming from the wavefunction of the nanoparticle. Basically, if we can improve the fit from this, by adding contact coupling, then we have some evidence that the wavefunction does extend off the nanoparticle.

### Lineshapes

If the hetergenous broadening is the dominant contributor (i.e., disorder in environments, movement during the experiment, etc) the bands will be Gaussian. Then we shape will be:

xxxx

Where we will need to take a guess at the width.

If the lifetime (i.e., $T_m$) dominates the broadening, then the shape will be Lorentzian.

xxx

It feels like we should be able to just straight calculate what this should be.  However, it might be convienent to include a scaling factor, to see how bad this Guess is.

#### Strain

For heterogeneous broadening, if we think that motion is important, then the motion during the experiment will have a greater effect near the particle than far away. This means that broadening will be larger when the coupling is larger (and depends more stronlgy on diestance) and so we can account for this by scaling the breadth by the coupling constant.

For homogeneous broadening... not totally sure.  I guess the question is if the lifetime depends on coupling at all.  I suppose stronger coupling means faster relaxation? But not sure.  Will need to talk with Alexey about this.

### Contact only

This will not have dipolar coupling. This is expected to be bad, since we should always have some dipolar coupling.

### Contact + Dipolar coupling

### Contact + Dipolar + strain broadening

## Size effects

Rather than broadening due to strain, it is possible that the size differences in the nanoparticles are responsible. This would have an effect on both contact and dipolar. So, if we account for the overall size distributions of our particles, we might get the broadening out right away.

### Contact + Size

### Dipolar + Size

### Dipolar + Contact + Size

### Dipolar + Contact + Size + Strain

## Results

Below are tables containing AIC values from fits.


### High density ligands

This is from the molecular dynamics simulations that have higher densities of ligands.


#### Mean size, no strain


| Lineshape  | dipolar only | contact only | both |
| ------------ | -------------- | -------------- | ------ |
| Lorentzian |              |              |      |
| Gaussian   |              |              |      |


#### Mean size, strain



| Lineshape  | dipolar only | contact only | both |
| ------------ | -------------- | -------------- | ------ |
| Lorentzian |              |              |      |
| Gaussian   |              |              |      |

#### Size distributions, no strain



| Lineshape  | dipolar only | contact only | both |
| ------------ | -------------- | -------------- | ------ |
| Lorentzian |              |              |      |
| Gaussian   |              |              |      |

#### Size distribution, strain



| Lineshape  | dipolar only | contact only | both |
| ------------ | -------------- | -------------- | ------ |
| Lorentzian |              |              |      |
| Gaussian   |              |              |      |
